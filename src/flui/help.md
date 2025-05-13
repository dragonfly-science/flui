
# Flui help

This application analyses FastQ files from a Nanopore sequencer.
It produces k-mer distributions from the FastQ reads to determine the subtype of the avian influenza virus present in each sample.
If the k-mer distributions match closely, the app will assign a subtype for the HA and NA segments.
Whether this subtype gets assigned depends on three thresholds:

* There must be a minimum number of k-mers found (to avoid premature matching).
* The matching score to particular subtype must be sufficient (to ensure the match is close).
* The matching score must be sufficiently bigger than the next biggest match (to ensure the match is not ambiguous).

## Locating the FastQ files

There are two ways that the app discovers FastQ files to process. When you start the application, you provide it with a parent folder and it will process:

1. **Pre-existing FastQ Files**: Any FASTQ files in sub-folders (regardless of the folder level) that match the Nanopore naming conventions. These files will be processed in random order.
2. **Incoming FastQ Files**: While the application is running,
  the app will monitor any subfolders for new FASTQ files that are placed there (by the Nanopore software as it processes the reads). These files will be processed in the order they arrive.

Each time a new FASTQ file for a particular barcode is processed, the scores and reads for that barcode are updated. This update happens in the background, so the effect may not be immediately evident.

## User interface

The interface contains five panes:

* Main pane: Each row has a run/barcode on it, and the columns show the status of that particular barcode. When the thresholds are met, the subtype is indicated. Until then, a question mark(`?`) is shown. If the subtype is ambiguous, then there will be two question marks (`??`). The other panes give more detail on why this ambiguity has occurred.
* NA and HA windows (on the right): The windows show the detail associated with the currently highlighted row in the main panel. Red writing shows which thresholds have not been met; green writing shows which thresholds have been passed.
* Events window (bottom): This window shows any new FASTQ files that are found, or processed, and which barcodes are being updated.
* Settings window (bottom right): This window shows the thresholds that are currently set (from the settings file, or the command line, or the defaults).

### Interactions

* Use the tab key and arrow keys to move around the interface.
* Sort the columns by clicking on the column header, or by using the shortcut keys shown at the bottom of the screen.

### Colour themes and screenshots (Ctrl-P)

The colour theme can be changed using the "Palette menu". Press Ctrl-P, and choose a theme that you like.
If you want to permanently select a theme, then you can set it in the `flui.toml` settings file.

You can also take a screenshot from the palette menu. Make sure to note down where the screenshot is saved. 
