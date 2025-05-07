
# Flui Help

This application analyses FastQ files from a Nanopore sequencer.
It produces kmer distributions from the FastQ reads to determine the subtype of the Avian Flu virus present in each sample.
If the kmer distributions match closely enough, the app will assign a subtype for the HA and NA segments.
Whether this gets assigned depends on three thresholds:

* There must be a minimum number of kmers found (to avoid premature matching).
* The matching score to particular subtype must be high enough (to ensure the match is good).
* The matching score must be sufficiently bigger than the next biggest match (to ensure the match is not ambiguous).

## Locating the FastQ Files

There are two ways that the app discovers FastQ files to process. When you start the application, you provide it with a parent folder and it will process:

1. **Pre-existing FastQ Files**: Any FASTQ files in sub-folders (however deep) that match the Nanopore naming conventions. These will be processed in random order.
2. **Incoming FastQ Files**: While the application is running,
  the app will monitor any subfolders for new FASTQ files that are placed there (by the Nanopore software as it processes the reads). These will be processed in the order they arrive.

Each time a new FASTQ file for a particular barcode is processed, the scores and reads for that barcode are updated. This happens in the background, so the effect may not be immediately seen.

## User Interface

The interface contains five panes:

* Main Pane: Each row has run/barcode on it, and the columns show the status of that particular barcode. When the thresholds are met, the subtype is indicated. Until then, we show a question mark(`?`). If the subtype is ambiguous, then there will be two question marks (`??`). The other panes give more detail on why this has occurred.
* NA and HA Windows (on the right): This shows the detail associated with the currently highlighted row in the main panel. Red writing shows which thresholds have not been met; green writing shows which thresholds are passed.
* Events Window (bottom): This shows any new FASTQ files that are found, or processed, and which barcodes are being updated.
* Settings Window (bottom right): This shows the thresholds that are currently set (from the settings file, or the command-line, or the defaults).

### Interactions

* Use the tab key and arrow keys to move around the interface.
* Sort the columns by clicking on the column header, or by using the shortcut keys shown at the bottom of the screen.

### Color Themes and Screenshots (Ctrl-P)

The color theme can be changed using the "Palette menu". Press Ctrl-P, and choose a theme that you like.
If you want to permanently set that theme, then you can do so in the `flui.toml` settings file.

You can also take a Screenshot from the palette menu. Make sure to note down where it puts the file.
