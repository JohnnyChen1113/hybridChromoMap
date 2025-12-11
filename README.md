# HybridChromoMap

A tool for drawing chromosome ancestry painting diagrams for hybrid species.

**Version: 0.1.0**

## Features

- Support for arbitrary species (user-defined chromosomes)
- Multiple origin coloring (not limited to 2 species, supports introgression)
- Variable ploidy support (haploid, diploid, aneuploid)
- Publication-quality output (SVG/PNG/PDF)
- Auto-generated colors when no color file is provided

## Installation

### Dependencies

- Python 3.8+
- matplotlib
- click

```bash
pip install matplotlib click
```

## Usage

```bash
# Basic usage (auto-generate colors)
python hybridchromomap.py -k karyotype.tsv -s segments.tsv -o output.png

# With custom colors
python hybridchromomap.py -k karyotype.tsv -s segments.tsv -c colors.tsv -o output.png

# Full options
python hybridchromomap.py \
    -k karyotype.tsv \
    -s segments.tsv \
    -c colors.tsv \
    -o output.svg \
    --sort length \
    --legend bottom \
    --width 14 \
    --dpi 300
```

### Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--karyotype` | `-k` | Karyotype definition file (TSV) | Required |
| `--segments` | `-s` | Segment origins file (TSV) | Required |
| `--colors` | `-c` | Origin colors file (TSV) | Auto-generate |
| `--out` | `-o` | Output file path | out.png |
| `--sort` | | Sort order (none/name/length) | none |
| `--legend` | | Legend position (right/bottom/none) | right |
| `--no-scale` | | Hide scale bar | False |
| `--width` | | Figure width in inches | 12 |
| `--dpi` | | PNG output resolution | 300 |

## Input File Formats

### Karyotype File (`karyotype.tsv`)

```tsv
#chrom	length	copy
chrI	500000	1
chrI	480000	2
chrII	800000	1
chrIII	600000	1
chrIII	590000	2
chrIII	595000	3
```

| Field | Type | Description |
|-------|------|-------------|
| chrom | string | Chromosome name |
| length | int | Length in bp |
| copy | int | Copy number (1, 2, 3...) |

### Segments File (`segments.tsv`)

```tsv
#chrom	copy	start	end	origin
chrI	1	0	250000	S_cerevisiae
chrI	1	250000	500000	S_paradoxus
chrI	2	0	200000	S_paradoxus
chrI	2	200000	480000	S_cerevisiae
```

| Field | Type | Description |
|-------|------|-------------|
| chrom | string | Chromosome name |
| copy | int | Copy number |
| start | int | Segment start position (bp) |
| end | int | Segment end position (bp) |
| origin | string | Origin identifier |

### Colors File (`colors.tsv`) - Optional

```tsv
#origin	color	label
S_cerevisiae	#E64B35	S. cerevisiae
S_paradoxus	#4DBBD5	S. paradoxus
unknown	#808080	Unknown
```

| Field | Type | Description |
|-------|------|-------------|
| origin | string | Origin identifier (matches segments file) |
| color | string | Hex color (e.g., #E64B35) |
| label | string | Display label for legend |

**Note:** If no colors file is provided, colors are auto-generated based on unique origins in the segments file.

## Examples

### Example 1: Basic Usage

A simple example with 4 chromosomes showing ancestry from multiple yeast species.

```bash
cd examples/example1_basic
python ../../hybridchromomap.py -k karyotype.tsv -s segments.tsv -c colors.tsv -o output.png
```

![Example 1 Output](examples/example1_basic/output.png)

### Example 2: Aneuploidy

Demonstrates variable ploidy: diploid (chrI, chrIV), haploid (chrII), and triploid (chrIII).

```bash
cd examples/example2_aneuploidy
python ../../hybridchromomap.py -k karyotype.tsv -s segments.tsv -o output.png
```

![Example 2 Output](examples/example2_aneuploidy/output.png)

### Example 3: S. bayanus Genome

Full 16-chromosome yeast genome with simulated ancestry painting.

```bash
cd examples/example3_sbay
python ../../hybridchromomap.py -k karyotype.tsv -s segments.tsv -o output.png --width 14
```

![Example 3 Output](examples/example3_sbay/output.png)

## Output Formats

Output format is determined by file extension:
- `.svg` - Vector graphics (editable)
- `.png` - Raster graphics (high resolution)
- `.pdf` - Vector graphics (publication-ready)

## License

MIT License

## Acknowledgments

Inspired by chromosome painting visualizations in hybrid species genomics research.
