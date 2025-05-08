# Flui 🦆🦠🧬

Flui is a command-line tool for **rapidly sub-typing avian influenza** viruses without doing full assemblies.

* Flui **identifies the HA and NA segments** of a virus from [Nanopore][nanopore] [FASTQ][fastq] files using kmer-based methods.
* Flui works with existing FASTQ files, but it can also monitor a folder for incoming FASTQ files, **providing real-time updates for an ongoing sequence run**.
* Flui runs in a terminal on any common platform, from your **Windows laptop to an SSH shell** on an HPC cluster.
* Flui has an **interactive user-interface** (TUI), showing continuous progress of the analysis.
* Flui uses a **simple, but robust metric**, to assign the subtype of the virus.

![A session of the interactive Flui interface](tui.png)

## Installation

> :warning: If you want to *develop*, rather than simply run `flui`, then see the section below for additional installation instructions.

### Recommended -- Install via UV

`flui` is a python package.
[UV][uv] is a modern tool that simplifies the installation of python packages, and is the recommended way to install `flui`.

To install `uv`, please follow the instructions from the [UV website][uv-install].

Once you gave `uv` installed, you can install `flui` with a single command:

```sh
uv tool install flui
```

You should now be able to run `flui --help`

### Alternative installation methods

You can install `flui` using any other traditional python methods (such as `pip`)
If you don't want to install directly from the internet, you can also install `flui` from a zip file.
These zip files are available from the [releases page](https://github.com/dragonfly-science/flui/releases).

## Usage

To get help on how to launch the UI, type `flui --help` into the terminal and press enter.
This will show the options available to you.
The two key things to provide are:

* `--run`: The path to a parent directory of the FastQ files.
  This should contain one or more runs, each containing multiple barcode sub-folders.
* `--ref`: This contains the reference genomes of HA and NA segments from different subtypes.

Typically, you will want to run a command like this:

```sh
flui --ref ref.fasta --run /path/to/fastq/files
```

## Test-driving for the impatient

Flui requires FastQ files from a Nanopore run, and a reference FASTA file.
If you don't have access to both of these, then you can try Flui out this way:

* Download the sample reference file [here][sample_ref] (created using the [NCBI virus data][ncbi]).
* Download sample FastQ files for Avian Flu from [this paper][sample_fastq].
  See the attached files section, and choose one or more of the zip files.
  Unzip the FastQ files into folders.

Once you have downloaded these you should unzip the FastQ downloads into a folder and then:

```sh
flui --ref reference-ncbi.fasta --run /folder/with/fastq
```

After a few moments, you should see the application start up and begin processing any existing FastQ files.

## Navigating the Application

Once you have started the application, you can navigate using the arrow keys and tab keys.
Detailed help is available inside the `flui` application.
Simply press the “h” button after starting the application.
You can also read it here: [help](src/flui/help.md).

## Saving Results

Flui saves the results of the analysis to both a CSV and a JSON file when it exits.
It saves the current state of the analysis, including the scores and the reads.
Note that this is saved *even if the analysis is incomplete*.
The files are saving in the current folder with date and time suffixes to prevent overwriting any existing files.

> :point_right: to avoid saving these files, start Flui with the `--no-dump` option.

## Configuration

The Flui app has several settings that can be changed, either at startup, or in a settings file.
The settings file must be called `flui.toml` and stored in the working directory.
Here you can set the kmer sizes, and the number of workers, and some UI colour options.
See the GitHub repository for an [example file][config].
Some settings can also be set on the command line (use `flui --help` to see these).
The settings are shown in the UI on the bottom right.

## How does it work?

Here is brief overview of Flui produces the scores for automatic sub-typing.

1. The `--ref` argument given on the command-line points to a [FASTA][fasta] file.
   This FASTA file contains the multiple reference sequences for each of the different subtypes (H1N1, H5N2, etc.).
   These sequences have both the subtype and segment number or type in the sequence header (i.e., HA/H1N1).
   We only use the HA and NA segments for sub-typing (others are ignored).
2. When Flui starts, it reads the FASTA file and, for each segment/subtype combination, it generates a kmer distribution.
   We store these distributions in memory.
3. Flui then reads in any existing FastQ files in folder and, for each read, it produces a kmer distribution.
   These kmer distributions are per run/barcode (we get this information from the file name).
   As more reads come in, we update the distribution for that run/barcode.
4. For each barcode distribution, we compare it to our set of reference distributions, and measure the *Jensen-Shannon Distance* (JSD) to each reference’s distribution.
   (The JSD is the square root of the [Jensen-Shannon Divergence][shannon], and is a proper [distance measure][metric]).
   The more the kmer distributions resemble each other, the lower the JSD.
5. We transform the JSD, to make it easier to interpret.
   First, we normalise it by dividing by the average JSD between all reference distributions.
   Call this the JSD*N*.
   Good matches will have JSDN values that fall below 1.0 (i.e. they are smaller than the distances between the references).
   To make this easier to interpret, we then take the complement of this value and multiply by 100: Matching Score = (1 - JSDN) \* 100.
6. The scores given in the UI are thus a *percentage reduction from expected kmer distribution distance*.
   Bigger values are better.
   Empirical tests show values of around 6.0 and above as typical for a good match.

## Development

For development, you’ll need to install the following dependencies:

* [uv][uv]
* [just][just]

Once `uv` is installed, you can use it to install some additional tools:

```sh
uv tool install ruff
uv tool install pyright
```

The `just` tool is used to run development-related tasks.
The `check` command runs all linting and tests, for example:

```sh
just check
```

At this stage, generating releases is not automated.

## Authorship and Funding

This project was developed by [Brett Calcott][brett] from [Dragonfly Data Science][dfly]
and [Ruy Jauregui][ruy] from [Biosecurity New Zealand][mpi].

Ruy managed reference development, sequence analysis, testing, and provided all bioinformatic guidance.
Brett was responsible for algorithm design and coding.

The project was funded by the [Biosecurity New Zealand
-- Tiakitanga Putaiao Aotearoa][mpi].

## License and Copyright

License [Apache 2.0][apache]

Copyright (c) 2024–2025 Dragonfly Data Science, Wellington, New Zealand.

[nanopore]: https://nanoporetech.com/platform/technology
[ncbi]: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide
[sample_ref]: https://github.com/dragonfly-science/flui/blob/main/sample/reference-ncbi.fasta
[sample_fastq]: https://www.sciencebase.gov/catalog/item/638a4df0d34ed907bf7907ea
[uv]: https://docs.astral.sh/uv/
[uv-install]: https://docs.astral.sh/uv/getting-started/installation/
[just]: https://github.com/casey/just
[fastq]: https://en.wikipedia.org/wiki/FASTQ_format
[shannon]: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
[metric]: https://en.wikipedia.org/wiki/Metric_space
[fasta]: https://en.wikipedia.org/wiki/FASTA_format
[brett]: https://github.com/brettc
[ruy]: https://github.com/ruy-jauregui
[mpi]: https://www.mpi.govt.nz/biosecurity/
[dfly]: https://www.dragonfly.co.nz
[apache]: https://www.apache.org/licenses/LICENSE-2.0
[config]: https://github.com/dragonfly-science/flui/blob/main/flui.toml
