from pathlib import Path
from typing import Any, ClassVar

from loguru import logger
from rich.console import RenderableType
from rich.text import Text
from textual import events, on
from textual._two_way_dict import TwoWayDict
from textual.app import App, ComposeResult
from textual.binding import Binding, BindingType
from textual.containers import Grid, Vertical, VerticalScroll
from textual.screen import ModalScreen
from textual.widgets import Button, DataTable, Footer, Label, Markdown, RichLog, Static
from textual.widgets.data_table import CellDoesNotExist, RowDoesNotExist, RowKey

from flui.dna import SegmentType
from flui.settings import Settings
from flui.subtype import Barcode, BarcodeProcessor, BarcodeSet, BarcodeUpdateKind


def to_kilos(count: int) -> str:
    if count == 0:
        return "0"
    if count < 1000:
        return "<.1k"
    return f"{count / 1000:.1f}k"


def build_text(settings: Settings) -> RenderableType:
    text = [
        f"No. Workers: [bold]{settings.workers}[/bold]",
        f"Min. Kmers: [bold]{settings.minimum_kmers}[/bold]",
        f"Score Threshold: [bold]{settings.minimum_score}[/bold]",
        f"Gap Threshold: [bold]{settings.minimum_gap}[/bold]",
        f"HA Kmer Size: [bold]{settings.ha_kmer_size}[/bold]",
        f"NA Kmer Size: [bold]{settings.na_kmer_size}[/bold]",
    ]
    return Text.from_markup("\n".join(text), justify="left")


class SettingsView(Static):
    def __init__(self, settings: Settings):
        super().__init__(build_text(settings), id="settings", classes="box")
        self.border_title = "Settings"


class EventView(RichLog):
    def __init__(self):
        super().__init__(id="events", classes="box", markup=True, max_lines=500)
        self.border_title = "Events"


class SegmentView(DataTable):
    def __init__(
        self,
        segment_type: SegmentType,
    ):
        super().__init__(id=segment_type.name, cursor_type="row", classes="box")
        self.cursor_type = "none"
        self.segment_type = segment_type
        self.border_title = segment_type.name
        self.can_focus = False
        self.build_columns()

    def build_columns(self):
        self.add_column(Text("Type", justify="center"), width=5, key="subtype")
        self.add_column(Text("Score", justify="center"), width=5, key="score")
        self.add_column(Text("Gap", justify="center"), width=5, key="gap")

    async def update(self, barcode: Barcode, settings: Settings):
        """Called when the selection of the barcode has changed."""
        self.clear()
        self.border_title = self.segment_type.name
        self.border_subtitle = barcode.key
        self.styles.border = ("round", "white")

        segment_match = (
            barcode.ha_match
            if self.segment_type is SegmentType.HA
            else barcode.na_match
        )
        if segment_match is None:
            return

        sc, fc = settings.success_color, settings.failure_color

        # We have a match, so update the border title and content.
        seg_count = (
            barcode.ha_kmer_count
            if self.segment_type is SegmentType.HA
            else barcode.na_kmer_count
        )
        kmers_in_kilos = to_kilos(seg_count)
        kmer_color = sc if seg_count >= settings.minimum_kmers else fc
        self.border_title = Text.from_markup(
            rf"{self.segment_type.name} [{kmer_color}]({kmers_in_kilos} kmers)"
        )

        last = None
        for i, (subtype, dist) in enumerate(segment_match.as_npdm()):
            if last is None:
                gap_text = ""
            else:
                gap = last - dist
                if i == 1:
                    gap_style = sc if gap >= settings.minimum_gap else fc
                else:
                    # Not interested in the gap after the first one.
                    gap_style = "dim"

                gap_text = Text(f"{gap:.3f}", justify="right", style=gap_style)

            if dist >= settings.minimum_score:
                style = sc
            elif i == 0:
                style = fc
            else:
                style = "dim"

            last = dist
            self.add_row(
                Text(subtype, justify="center", style=style),
                Text(f"{dist:.3f}", justify="right", style=style),
                gap_text,
                key=subtype,
            )

        if "?" not in barcode.get_assigned(self.segment_type, settings):
            self.styles.border = ("heavy", sc)


class BarcodeView(DataTable):
    def __init__(self, barcodes: BarcodeSet, settings: Settings):
        super().__init__(id="barcode", cursor_type="row", classes="box left")
        self.border_title = "Barcodes"
        self.barcodes = barcodes
        self.settings = settings
        self.build_columns()
        self.header_height = 2
        self.sort_type = "default"
        self.sort_reverse = False

    @on(DataTable.HeaderSelected)
    async def header_selected(self, selected: DataTable.HeaderSelected):
        if selected.column_key.value is None:
            return
        await self.change_sort(selected.column_key.value)

    async def change_sort(self, column_key: str):
        if column_key in {"run", "flowcell", "barcode"}:
            sort_type = "default"
        else:
            sort_type = column_key

        if self.sort_type == sort_type:
            self.sort_reverse = not self.sort_reverse

        self.sort_type = sort_type
        self.sort_rows()

    def sort_rows(self):
        selected = self.get_selected_key()

        def key_wrapper(row: tuple[RowKey, Any]) -> Any:
            row_key, _ = row
            assert row_key.value is not None
            barcode = self.barcodes.barcodes[row_key.value]
            default = (barcode.run_id, barcode.flow_cell_id, barcode.barcode_id)
            # subtype = barcode.self.get_subtype(barcode)
            sorter = {
                "files": barcode.file_count,
                "done": barcode.done_count,
                "ha": barcode.get_assigned(SegmentType.HA, self.settings),
                "na": barcode.get_assigned(SegmentType.NA, self.settings),
                "reads": barcode.reads_count,
                "ha_kmers": barcode.ha_kmer_count,
                "na_kmers": barcode.na_kmer_count,
                "ha_score": barcode.ha_score,
                "na_score": barcode.na_score,
            }
            return sorter.get(self.sort_type, default)

        # NOTE: This is a hack taken from the internals of DataTable
        # As we can't sort the rows directly by row key. Which is what we want.
        # TODO: Submit a PR?
        ordered_rows = sorted(
            self._data.items(), key=key_wrapper, reverse=self.sort_reverse
        )
        self._row_locations = TwoWayDict(
            {row_key: new_index for new_index, (row_key, _) in enumerate(ordered_rows)}
        )
        self._update_count += 1
        self.refresh()
        self.set_selected_key(selected)

    def get_selected_key(self) -> str | None:
        try:
            cell_key = self.coordinate_to_cell_key(self.cursor_coordinate)
        except CellDoesNotExist:
            return None
        return cell_key.row_key.value

    def set_selected_key(self, key: str | None):
        if key is None:
            return

        try:
            row = self.get_row_index(key)
        except RowDoesNotExist:
            pass
        else:
            self.move_cursor(row=row)

    def build_columns(self):
        self.add_column(Text("\nRun Id", justify="center"), width=9, key="runid")
        self.add_column(Text("\nFlowcell", justify="center"), width=9, key="flowcell")
        self.add_column(Text("Bar-\ncode", justify="center"), width=5, key="barcode")
        self.add_column(Text("\nHA", justify="center"), width=3, key="ha")
        self.add_column(Text("\nNA", justify="center"), width=3, key="na")
        self.add_column(Text("FASTQ\nTotal", justify="center"), width=5, key="files")
        self.add_column(Text("FASTQ\nDone", justify="center"), width=5, key="done")
        self.add_column(Text("\nReads", justify="center"), width=7, key="reads")
        self.add_column(Text("HA\nKmers", justify="center"), width=6, key="ha_kmers")
        self.add_column(Text("NA\nKmers", justify="center"), width=6, key="na_kmers")
        self.add_column(Text("HA\nScore", justify="center"), width=5, key="ha_score")
        self.add_column(Text("NA\nScore", justify="center"), width=5, key="na_score")

    def build_record(self, barcode: Barcode) -> dict[str, Text]:
        sc = self.settings.success_color

        def style_subtype(segment_type: SegmentType) -> Text:
            st = barcode.get_assigned(segment_type, self.settings)
            style = "" if "?" in st else sc
            return Text(st, justify="center", style=style)

        def style_score(score: float) -> Text:
            style = sc if score >= self.settings.minimum_score else ""
            return Text(f"{score:.2f}", justify="right", style=style)

        def style_count(count: int) -> Text:
            style = sc if count >= self.settings.minimum_kmers else ""
            return Text(f"{to_kilos(count)}", justify="right", style=style)

        return dict(
            runid=Text(barcode.run_id, style="bold", justify="left"),
            flowcell=Text(barcode.flow_cell_id, style="bold", justify="left"),
            barcode=Text(barcode.barcode_id, style="bold", justify="center"),
            ha=style_subtype(SegmentType.HA),
            na=style_subtype(SegmentType.NA),
            files=Text(f"{barcode.file_count}", justify="center"),
            done=Text(f"{barcode.done_count}", justify="center"),
            reads=Text(f"{to_kilos(barcode.reads_count)}", justify="right"),
            ha_kmers=style_count(barcode.ha_kmer_count),
            na_kmers=style_count(barcode.na_kmer_count),
            ha_score=style_score(barcode.ha_score),
            na_score=style_score(barcode.na_score),
        )

    async def load(self):
        """Initial loading."""
        for barcode in self.barcodes.barcodes.values():
            columns = self.build_record(barcode)
            self.add_row(*columns.values(), key=barcode.key)

        self.sort_rows()

    async def update_barcode(self, barcode: Barcode):
        # Have we already got this barcode.
        try:
            self.get_row(barcode.key)
            has_row = True
        except RowDoesNotExist:
            has_row = False

        columns = self.build_record(barcode)

        if not has_row:
            # Add a new row. Record where we are... (using key)
            key = self.get_selected_key()
            self.add_row(*columns.values(), key=barcode.key)
            # Reset the same key
            self.set_selected_key(key)
        else:
            # Update existing row.
            for k, text in columns.items():
                self.update_cell(barcode.key, k, text)

        self.sort_rows()


class QuitScreen(ModalScreen[bool]):
    """Screen to run."""

    BINDINGS: ClassVar[list[BindingType]] = [
        Binding("y", "run", "Yes", show=False),
        Binding("escape,n", "cancel", "No", show=False),
    ]

    def compose(self) -> ComposeResult:
        yield Grid(
            Label("Do you really want to quit?", id="question"),
            Button("[bold underline]Y[/bold underline]es", variant="error", id="yep"),
            Button("[bold underline]N[/bold underline]o", variant="primary", id="nope"),
            id="dialog",
        )

    async def on_mount(self):
        button = self.query_one("#nope", Button)
        button.focus()

    async def action_run(self):
        self.dismiss(True)

    async def action_cancel(self):
        self.dismiss(False)

    @on(Button.Pressed, "#yep")
    def handle_run(self):
        self.dismiss(True)

    @on(Button.Pressed, "#nope")
    def handle_cancel(self):
        self.dismiss(False)


# Gives us some pretty arrows that mean "sort"
UP_DOWN_CHARS = "\u2191\u2193"


class FluiApp(App):
    TITLE = "Flui Subtyping"
    CSS_PATH = "app.scss"
    BINDINGS = [  # noqa: RUF012
        ("h", "help_screen", "Help"),
        ("ctrl+q", "quit", "Quit"),
        ("b", "sort_by('default')", UP_DOWN_CHARS + "Run/Barcode"),
        ("s", "sort_by('ha_score')", UP_DOWN_CHARS + "HA"),
        ("S", "sort_by('na_score')", UP_DOWN_CHARS + "NA"),
        ("r", "sort_by('reads')", UP_DOWN_CHARS + "Reads"),
        ("k", "sort_by('ha_kmers')", UP_DOWN_CHARS + "HA-Kmers"),
        ("K", "sort_by('na_kmers')", UP_DOWN_CHARS + "NA-Kmers"),
    ]

    def __init__(self, barcode_set: BarcodeSet, settings: Settings):
        super().__init__()
        self.sub_title = str(barcode_set.root)
        self.barcode_set = barcode_set
        self.processor = BarcodeProcessor(
            barcode_set=self.barcode_set, max_workers=settings.workers
        )
        self.settings = settings
        self.quitting = False

    def on_ready(self):
        self.set_interval(0.5, self.process_updates)

    async def action_sort_by(self, sort_type: str):
        await self.barcode_view.change_sort(sort_type)

    async def action_help_screen(self) -> None:
        await self.push_screen(HelpScreen(id="help_screen"))

    def compose(self) -> ComposeResult:
        """Create child widgets for the app."""
        self.barcode_view = BarcodeView(self.barcode_set, self.settings)
        self.info_view = SettingsView(self.settings)
        self.event_view = EventView()
        self.ha_view = SegmentView(SegmentType.HA)
        self.na_view = SegmentView(SegmentType.NA)

        yield self.barcode_view
        yield self.ha_view
        yield self.na_view
        yield self.event_view
        yield self.info_view
        yield Footer()

    async def on_mount(self):
        self.theme = self.settings.theme
        await self.barcode_view.load()
        self.processor.start()

    async def action_quit(self):
        def check_submit(result: bool | None):
            self.quitting = False
            if result:
                # Shut down background processor.
                self.processor.shutdown()
                self.exit()

        # Prevent doubling up of the quit screen.
        # Not sure why this happens.
        if self.quitting:
            return

        self.quitting = True
        await self.push_screen(QuitScreen(), check_submit)

    async def process_updates(self):
        """Read in any of the completions from background processing."""
        updates = await self.processor.get_updates()
        if updates:
            logger.info(f"TUI processing {len(updates)} updates")
        for update in updates:
            logger.info(f"TUI processing update: {update}")
            await self.barcode_view.update_barcode(update.barcode)
            logger.info(f"Writing update to event view: {update}")
            self.event_view.write(update)

            if update.kind is BarcodeUpdateKind.UPDATE_MATCHES:
                bc = update.barcode
                if self.barcode_view.get_selected_key() == bc.key:
                    await self.ha_view.update(bc, self.settings)
                    await self.na_view.update(bc, self.settings)

    @on(DataTable.RowHighlighted, "#barcode")
    async def change_barcode(self):
        barcode_name = self.barcode_view.get_selected_key()
        if barcode_name is not None:
            bc = self.barcode_set.barcodes[barcode_name]
            await self.ha_view.update(bc, self.settings)
            await self.na_view.update(bc, self.settings)


class VerticalSuppressClicks(Vertical):
    def on_click(self, message: events.Click) -> None:
        message.stop()


class HelpScreen(ModalScreen):
    def compose(self) -> ComposeResult:
        markdown_path = Path(__file__).parent / "help.md"
        with markdown_path.open("r") as f:
            markdown = f.read()

        with VerticalSuppressClicks(id="help_outer"):
            with VerticalScroll(id="help_inner"):
                yield Markdown(markdown=markdown)
            yield Static(
                "Scroll with arrows. Press any other key to continue.", id="help_footer"
            )

    def on_mount(self) -> None:
        self.body = self.query_one("#help_inner")

    def on_key(self, event: events.Key) -> None:
        event.stop()
        if event.key == "up":
            self.body.scroll_up()
        elif event.key == "down":
            self.body.scroll_down()
        elif event.key == "left":
            self.body.scroll_left()
        elif event.key == "right":
            self.body.scroll_right()
        elif event.key == "pageup":
            self.body.scroll_page_up()
        elif event.key == "pagedown":
            self.body.scroll_page_down()
        else:
            self.app.pop_screen()

    def on_click(self) -> None:
        self.app.pop_screen()
