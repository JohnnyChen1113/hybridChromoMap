#!/usr/bin/env python3
"""
HybridChromoMap - Chromosome ancestry painting visualization tool

A tool for drawing hybrid species chromosome ancestry painting diagrams.
Supports arbitrary species, multiple origins (including introgression),
and variable ploidy levels.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import click
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.collections import PatchCollection


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class Origin:
    """Represents a genomic origin (species/source) with its display properties."""
    name: str
    color: str  # Hex color like #E63946
    label: str  # Display label for legend

    def to_rgb(self) -> Tuple[float, float, float]:
        """Convert hex color to RGB tuple (0-1 range)."""
        hex_color = self.color.lstrip('#')
        return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


@dataclass
class Segment:
    """Represents a genomic segment with a specific origin."""
    chrom: str
    copy: int
    start: int
    end: int
    origin: str  # References Origin.name

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class ChromosomeCopy:
    """Represents one copy of a chromosome."""
    chrom: str
    copy: int
    length: int
    segments: List[Segment] = field(default_factory=list)

    def add_segment(self, segment: Segment):
        self.segments.append(segment)

    def sort_segments(self):
        """Sort segments by start position."""
        self.segments.sort(key=lambda s: s.start)


@dataclass
class Chromosome:
    """Represents a chromosome with potentially multiple copies."""
    name: str
    copies: Dict[int, ChromosomeCopy] = field(default_factory=dict)

    def add_copy(self, copy_num: int, length: int):
        self.copies[copy_num] = ChromosomeCopy(
            chrom=self.name,
            copy=copy_num,
            length=length
        )

    @property
    def max_length(self) -> int:
        """Return the maximum length among all copies."""
        return max(c.length for c in self.copies.values()) if self.copies else 0

    @property
    def ploidy(self) -> int:
        """Return the number of copies."""
        return len(self.copies)


@dataclass
class Karyotype:
    """Manages all chromosome information."""
    chromosomes: Dict[str, Chromosome] = field(default_factory=dict)
    chrom_order: List[str] = field(default_factory=list)  # Preserve input order

    def add_chromosome_copy(self, chrom_name: str, copy_num: int, length: int):
        if chrom_name not in self.chromosomes:
            self.chromosomes[chrom_name] = Chromosome(name=chrom_name)
            self.chrom_order.append(chrom_name)
        self.chromosomes[chrom_name].add_copy(copy_num, length)

    def get_copy(self, chrom_name: str, copy_num: int) -> Optional[ChromosomeCopy]:
        if chrom_name in self.chromosomes:
            return self.chromosomes[chrom_name].copies.get(copy_num)
        return None

    @property
    def max_length(self) -> int:
        """Return the maximum chromosome length across all chromosomes."""
        return max(c.max_length for c in self.chromosomes.values()) if self.chromosomes else 0

    def get_ordered_chromosomes(self, sort_by: str = 'none') -> List[str]:
        """Return chromosome names in specified order."""
        if sort_by == 'name':
            return sorted(self.chrom_order)
        elif sort_by == 'length':
            return sorted(self.chrom_order,
                         key=lambda n: self.chromosomes[n].max_length,
                         reverse=True)
        else:  # 'none' - preserve input order
            return self.chrom_order


# =============================================================================
# File Parsers
# =============================================================================

def parse_karyotype(filepath: Path) -> Karyotype:
    """
    Parse karyotype TSV file.
    Format: #chrom  length  copy
    """
    karyotype = Karyotype()

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                raise ValueError(
                    f"Karyotype file line {line_num}: expected 3 columns "
                    f"(chrom, length, copy), got {len(parts)}"
                )

            chrom = parts[0]
            try:
                length = int(parts[1])
                copy = int(parts[2])
            except ValueError as e:
                raise ValueError(
                    f"Karyotype file line {line_num}: invalid integer value: {e}"
                )

            karyotype.add_chromosome_copy(chrom, copy, length)

    return karyotype


def parse_segments(filepath: Path, karyotype: Karyotype) -> None:
    """
    Parse segments TSV file and add segments to karyotype.
    Format: #chrom  copy  start  end  origin
    """
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 5:
                raise ValueError(
                    f"Segments file line {line_num}: expected 5 columns "
                    f"(chrom, copy, start, end, origin), got {len(parts)}"
                )

            chrom = parts[0]
            try:
                copy = int(parts[1])
                start = int(parts[2])
                end = int(parts[3])
            except ValueError as e:
                raise ValueError(
                    f"Segments file line {line_num}: invalid integer value: {e}"
                )
            origin = parts[4]

            segment = Segment(chrom=chrom, copy=copy, start=start, end=end, origin=origin)

            chrom_copy = karyotype.get_copy(chrom, copy)
            if chrom_copy is None:
                raise ValueError(
                    f"Segments file line {line_num}: chromosome '{chrom}' copy {copy} "
                    f"not found in karyotype"
                )
            chrom_copy.add_segment(segment)

    # Sort segments for each copy
    for chrom in karyotype.chromosomes.values():
        for copy in chrom.copies.values():
            copy.sort_segments()


def parse_color(color_str: str) -> str:
    """
    Parse color string and return hex format (#RRGGBB).

    Supported formats:
    - Hex: #E64B35, #e64b35
    - RGB tuple: (230, 75, 53), (230,75,53), 230,75,53
    - Color names: red, blue, green, etc. (matplotlib named colors)
    """
    color_str = color_str.strip()

    # Format 1: Hex color (#RRGGBB or #RGB)
    if color_str.startswith('#'):
        hex_color = color_str[1:]
        if len(hex_color) == 6:
            # Validate hex characters
            try:
                int(hex_color, 16)
                return color_str.upper()
            except ValueError:
                pass
        elif len(hex_color) == 3:
            # Short hex format #RGB -> #RRGGBB
            try:
                int(hex_color, 16)
                return f"#{hex_color[0]*2}{hex_color[1]*2}{hex_color[2]*2}".upper()
            except ValueError:
                pass

    # Format 2: RGB tuple - (R, G, B) or R,G,B
    # Remove parentheses and split
    rgb_str = color_str.strip('()').replace(' ', '')
    if ',' in rgb_str:
        parts = rgb_str.split(',')
        if len(parts) == 3:
            try:
                r, g, b = [int(p) for p in parts]
                if all(0 <= v <= 255 for v in (r, g, b)):
                    return f"#{r:02X}{g:02X}{b:02X}"
            except ValueError:
                pass

    # Format 3: Named color (matplotlib)
    try:
        rgb = mcolors.to_rgb(color_str)
        r, g, b = [int(v * 255) for v in rgb]
        return f"#{r:02X}{g:02X}{b:02X}"
    except ValueError:
        pass

    raise ValueError(
        f"Invalid color format: '{color_str}'. "
        f"Supported formats: #RRGGBB, (R,G,B), or color name (red, blue, etc.)"
    )


def parse_origins(filepath: Path) -> Dict[str, Origin]:
    """
    Parse origins TSV file.
    Format: #origin  color  label

    Color can be specified as:
    - Hex: #E64B35
    - RGB: (230, 75, 53) or 230,75,53
    - Name: red, blue, green, etc.
    """
    origins = {}

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                raise ValueError(
                    f"Origins file line {line_num}: expected 3 columns "
                    f"(origin, color, label), got {len(parts)}"
                )

            name = parts[0]
            color_str = parts[1]
            label = parts[2]

            try:
                color = parse_color(color_str)
            except ValueError as e:
                raise ValueError(f"Origins file line {line_num}: {e}")

            origins[name] = Origin(name=name, color=color, label=label)

    # Ensure 'unknown' origin exists
    if 'unknown' not in origins:
        origins['unknown'] = Origin(name='unknown', color='#808080', label='Unknown')

    return origins


def validate_data(karyotype: Karyotype, origins: Dict[str, Origin]) -> List[str]:
    """Validate that all segment origins exist in origins dict."""
    errors = []

    for chrom in karyotype.chromosomes.values():
        for copy in chrom.copies.values():
            for segment in copy.segments:
                if segment.origin not in origins:
                    errors.append(
                        f"Segment origin '{segment.origin}' for {segment.chrom} "
                        f"copy {segment.copy} not found in origins file"
                    )

    return errors


def generate_origin_colors(karyotype: Karyotype) -> Dict[str, Origin]:
    """
    Generate automatic colors for each unique origin in the segments.
    Extracts all unique origin values from segments and assigns distinct colors.
    """
    # NPG (Nature Publishing Group) inspired color palette
    color_palette = [
        '#E64B35',  # Red (主色)
        '#4DBBD5',  # Blue (蓝色)
        '#00A087',  # Teal/黄绿色
        '#3C5488',  # Dark Blue/绿色
        '#F39B7F',  # Peach/黄橙色
        '#8491B4',  # Lavender/紫色
        '#91D1C2',  # Mint/橙色
        '#DC0000',  # Crimson/粉红色
        '#7E6148',  # Brown/浅蓝色
        '#000000',  # Black/黑色
    ]

    # Collect all unique origins from segments
    unique_origins = set()
    for chrom in karyotype.chromosomes.values():
        for copy in chrom.copies.values():
            for segment in copy.segments:
                unique_origins.add(segment.origin)

    # Remove 'unknown' from the set if present (we'll add it separately)
    unique_origins.discard('unknown')

    # Sort origins for consistent coloring
    sorted_origins = sorted(unique_origins)

    origins = {}
    for i, origin_name in enumerate(sorted_origins):
        color = color_palette[i % len(color_palette)]
        # Create a display label (replace underscores with spaces, add italic formatting hint)
        label = origin_name.replace('_', ' ')
        origins[origin_name] = Origin(
            name=origin_name,
            color=color,
            label=label
        )

    # Add unknown as fallback with gray color
    origins['unknown'] = Origin(name='unknown', color='#808080', label='Unknown')

    return origins


# =============================================================================
# Renderer
# =============================================================================

class ChromoMapRenderer:
    """Renders chromosome ancestry painting maps."""

    def __init__(
        self,
        karyotype: Karyotype,
        origins: Dict[str, Origin],
        fig_width: float = 12.0,
        chrom_height: float = 0.4,
        font_size: float = 10.0,
        same_chrom_spacing: float = 0.15,
        diff_chrom_spacing: float = 0.5,
        left_margin: float = 1.5,
        right_margin: float = 2.5,
        top_margin: float = 0.5,
        bottom_margin: float = 1.0,
    ):
        self.karyotype = karyotype
        self.origins = origins
        self.fig_width = fig_width
        self.chrom_height = chrom_height
        self.font_size = font_size
        self.same_chrom_spacing = same_chrom_spacing
        self.diff_chrom_spacing = diff_chrom_spacing
        self.left_margin = left_margin
        self.right_margin = right_margin
        self.top_margin = top_margin
        self.bottom_margin = bottom_margin

        # Calculate scale: bp to inches
        self.max_bp = karyotype.max_length
        self.plot_width = fig_width - left_margin - right_margin
        self.bp_to_inch = self.plot_width / self.max_bp if self.max_bp > 0 else 1

    def _calculate_fig_height(self, sort_by: str) -> float:
        """Calculate figure height based on number of chromosome copies."""
        total_height = self.top_margin + self.bottom_margin

        chrom_names = self.karyotype.get_ordered_chromosomes(sort_by)
        for i, chrom_name in enumerate(chrom_names):
            chrom = self.karyotype.chromosomes[chrom_name]
            num_copies = len(chrom.copies)

            # Add height for chromosome copies
            total_height += num_copies * self.chrom_height
            # Add spacing between copies of same chromosome
            total_height += (num_copies - 1) * self.same_chrom_spacing

            # Add spacing to next chromosome (except for last)
            if i < len(chrom_names) - 1:
                total_height += self.diff_chrom_spacing

        return total_height

    def _draw_capsule(
        self,
        ax: plt.Axes,
        x: float,
        y: float,
        width: float,
        height: float,
        facecolor: str = '#CCCCCC',
        edgecolor: str = 'black',
        linewidth: float = 0.5,
        zorder: int = 1
    ):
        """Draw a capsule (rounded rectangle with half-circle ends)."""
        # Use FancyBboxPatch with round corners
        rounding = min(height / 2, width / 2)
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle=f"round,pad=0,rounding_size={rounding}",
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder
        )
        ax.add_patch(patch)
        return patch

    def _draw_segment(
        self,
        ax: plt.Axes,
        x: float,
        y: float,
        width: float,
        height: float,
        color: str,
        is_start: bool = False,
        is_end: bool = False,
        chrom_length_inch: float = 0,
        zorder: int = 2
    ):
        """Draw a segment, with rounded ends if at chromosome boundaries."""
        rounding = height / 2

        if is_start and is_end:
            # Entire chromosome is one segment - full capsule
            patch = FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle=f"round,pad=0,rounding_size={rounding}",
                facecolor=color,
                edgecolor='none',
                zorder=zorder
            )
        elif is_start:
            # Left end rounded
            patch = FancyBboxPatch(
                (x, y),
                width + rounding,
                height,
                boxstyle=f"round,pad=0,rounding_size={rounding}",
                facecolor=color,
                edgecolor='none',
                zorder=zorder,
                clip_on=True
            )
            # Clip right side
            clip_rect = Rectangle((x, y), width, height)
            patch.set_clip_path(clip_rect, transform=ax.transData)
        elif is_end:
            # Right end rounded
            patch = FancyBboxPatch(
                (x - rounding, y),
                width + rounding,
                height,
                boxstyle=f"round,pad=0,rounding_size={rounding}",
                facecolor=color,
                edgecolor='none',
                zorder=zorder,
                clip_on=True
            )
            # Clip left side
            clip_rect = Rectangle((x, y), width, height)
            patch.set_clip_path(clip_rect, transform=ax.transData)
        else:
            # Middle segment - rectangular
            patch = Rectangle(
                (x, y),
                width,
                height,
                facecolor=color,
                edgecolor='none',
                zorder=zorder
            )

        ax.add_patch(patch)
        return patch

    def _draw_scale_bar(self, ax: plt.Axes, y_pos: float):
        """Draw a scale bar at the bottom of the figure."""
        # Determine appropriate scale unit
        if self.max_bp >= 1_000_000:
            scale_unit = 1_000_000
            unit_name = "Mb"
        else:
            scale_unit = 1_000
            unit_name = "kb"

        max_val = self.max_bp / scale_unit

        # Choose tick interval to have ~5-8 ticks maximum
        # Find a nice interval that gives reasonable number of ticks
        nice_intervals = [0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000]
        tick_interval = nice_intervals[0]
        for interval in nice_intervals:
            num_ticks = max_val / interval
            if 4 <= num_ticks <= 8:
                tick_interval = interval
                break
            elif num_ticks < 4:
                # Previous interval was better, but this one gives too few ticks
                break
            tick_interval = interval

        # Draw axis line
        x_start = self.left_margin
        x_end = self.left_margin + self.plot_width
        ax.plot([x_start, x_end], [y_pos, y_pos], 'k-', linewidth=1, zorder=10)

        # Draw ticks and labels
        tick_height = 0.1
        tick_val = 0
        while tick_val <= max_val + 0.001:  # Small tolerance for float comparison
            x = self.left_margin + (tick_val * scale_unit) * self.bp_to_inch
            if x <= x_end + 0.01:
                ax.plot([x, x], [y_pos, y_pos - tick_height], 'k-', linewidth=1, zorder=10)
                # Format label nicely
                if tick_val == int(tick_val):
                    label = f"{int(tick_val)}"
                else:
                    label = f"{tick_val:g}"
                ax.text(x, y_pos - tick_height - 0.1, label,
                       ha='center', va='top', fontsize=9, zorder=10)
            tick_val += tick_interval

        # Add unit label
        ax.text(x_end + 0.1, y_pos, unit_name, ha='left', va='center', fontsize=10, zorder=10)

    def _draw_legend(
        self,
        ax: plt.Axes,
        position: str,
        fig_height: float
    ):
        """Draw the legend."""
        if position == 'none':
            return

        # Get unique origins that are actually used
        used_origins = set()
        for chrom in self.karyotype.chromosomes.values():
            for copy in chrom.copies.values():
                for segment in copy.segments:
                    used_origins.add(segment.origin)

        # Create legend entries
        handles = []
        labels = []
        for origin_name in sorted(used_origins):
            if origin_name in self.origins:
                origin = self.origins[origin_name]
                handle = mpatches.Patch(facecolor=origin.color, edgecolor='none')
                handles.append(handle)
                labels.append(origin.label)

        if not handles:
            return

        if position == 'right':
            legend = ax.legend(
                handles, labels,
                loc='center left',
                bbox_to_anchor=(1.02, 0.5),
                frameon=True,
                fontsize=10,
                title='Origin'
            )
        elif position == 'bottom':
            legend = ax.legend(
                handles, labels,
                loc='upper center',
                bbox_to_anchor=(0.5, -0.15),
                ncol=min(len(handles), 4),
                frameon=True,
                fontsize=10,
                title='Origin'
            )

    def render(
        self,
        output_path: Path,
        sort_by: str = 'none',
        legend_position: str = 'right',
        show_scale: bool = True,
        dpi: int = 300
    ):
        """Render the chromosome map to a file."""
        fig_height = self._calculate_fig_height(sort_by)

        fig, ax = plt.subplots(figsize=(self.fig_width, fig_height))
        ax.set_xlim(0, self.fig_width)
        ax.set_ylim(0, fig_height)
        ax.axis('off')

        # Start from top
        current_y = fig_height - self.top_margin

        chrom_names = self.karyotype.get_ordered_chromosomes(sort_by)

        for i, chrom_name in enumerate(chrom_names):
            chrom = self.karyotype.chromosomes[chrom_name]

            # Sort copies by copy number
            sorted_copies = sorted(chrom.copies.keys())

            for j, copy_num in enumerate(sorted_copies):
                copy = chrom.copies[copy_num]

                # Calculate y position for this copy
                y = current_y - self.chrom_height

                # Calculate chromosome length in inches
                chrom_length_inch = copy.length * self.bp_to_inch

                # Draw chromosome outline (capsule)
                self._draw_capsule(
                    ax,
                    self.left_margin,
                    y,
                    chrom_length_inch,
                    self.chrom_height,
                    facecolor='#E0E0E0',
                    edgecolor='#404040',
                    linewidth=0.5,
                    zorder=1
                )

                # Draw segments
                rounding = self.chrom_height / 2

                for seg_idx, segment in enumerate(copy.segments):
                    seg_x = self.left_margin + segment.start * self.bp_to_inch
                    seg_width = segment.length * self.bp_to_inch

                    # Check if this is the start/end of chromosome
                    is_start = (segment.start == 0)
                    is_end = (segment.end >= copy.length - 1)

                    # Get segment color
                    color = self.origins.get(segment.origin, self.origins['unknown']).color

                    # For proper rounded ends, draw as capsule if at edges
                    if is_start and is_end:
                        # Single segment covers entire chromosome - full capsule
                        self._draw_capsule(
                            ax, seg_x, y, seg_width, self.chrom_height,
                            facecolor=color, edgecolor='none', zorder=2
                        )
                    elif is_start:
                        # First segment - round left edge, flat right edge
                        # Draw a full rounded box then cover right rounded part
                        patch = FancyBboxPatch(
                            (seg_x, y),
                            seg_width + rounding,
                            self.chrom_height,
                            boxstyle=f"round,pad=0,rounding_size={rounding}",
                            facecolor=color,
                            edgecolor='none',
                            zorder=2
                        )
                        ax.add_patch(patch)
                        # Cover the right rounded part
                        cover = Rectangle(
                            (seg_x + seg_width, y),
                            rounding,
                            self.chrom_height,
                            facecolor='#E0E0E0',
                            edgecolor='none',
                            zorder=1
                        )
                        ax.add_patch(cover)
                    elif is_end:
                        # Last segment - flat left edge, round right edge
                        # Draw rectangle first, then add a rounded cap on the right
                        # Main rectangle part
                        rect = Rectangle(
                            (seg_x, y),
                            seg_width - rounding,
                            self.chrom_height,
                            facecolor=color,
                            edgecolor='none',
                            zorder=2
                        )
                        ax.add_patch(rect)
                        # Right rounded cap (half circle)
                        from matplotlib.patches import Arc, Wedge, Circle
                        cap_center_x = seg_x + seg_width - rounding
                        cap_center_y = y + self.chrom_height / 2
                        # Use a wedge for the right semicircle
                        cap = Wedge(
                            (cap_center_x, cap_center_y),
                            rounding,
                            -90, 90,  # Right semicircle
                            facecolor=color,
                            edgecolor='none',
                            zorder=2
                        )
                        ax.add_patch(cap)
                    else:
                        # Middle segment - plain rectangle
                        patch = Rectangle(
                            (seg_x, y),
                            seg_width,
                            self.chrom_height,
                            facecolor=color,
                            edgecolor='none',
                            zorder=2
                        )
                        ax.add_patch(patch)

                # Draw label
                label = f"{chrom_name}-{copy_num}"
                ax.text(
                    self.left_margin - 0.1,
                    y + self.chrom_height / 2,
                    label,
                    ha='right',
                    va='center',
                    fontsize=self.font_size,
                    fontfamily='monospace'
                )

                # Update y position
                current_y = y
                if j < len(sorted_copies) - 1:
                    current_y -= self.same_chrom_spacing

            # Add larger spacing before next chromosome
            if i < len(chrom_names) - 1:
                current_y -= self.diff_chrom_spacing

        # Draw scale bar
        if show_scale:
            scale_y = self.bottom_margin * 0.6
            self._draw_scale_bar(ax, scale_y)

        # Draw legend
        self._draw_legend(ax, legend_position, fig_height)

        # Save figure
        plt.tight_layout()

        suffix = output_path.suffix.lower()
        if suffix == '.svg':
            plt.savefig(output_path, format='svg', bbox_inches='tight')
        elif suffix == '.pdf':
            plt.savefig(output_path, format='pdf', bbox_inches='tight')
        else:
            plt.savefig(output_path, format='png', dpi=dpi, bbox_inches='tight')

        plt.close(fig)


# =============================================================================
# CLI
# =============================================================================

@click.command()
@click.option(
    '-k', '--karyotype',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Karyotype definition file (TSV)'
)
@click.option(
    '-s', '--segments',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Segment origins file (TSV)'
)
@click.option(
    '-c', '--colors',
    required=False,
    default=None,
    type=click.Path(exists=True, path_type=Path),
    help='Origin colors definition file (TSV). If not provided, auto-color by chromosome.'
)
@click.option(
    '-o', '--out',
    default='out.png',
    type=click.Path(path_type=Path),
    help='Output file path (SVG/PNG/PDF)'
)
@click.option(
    '--sort',
    type=click.Choice(['none', 'name', 'length']),
    default='none',
    help='Chromosome sort order'
)
@click.option(
    '--legend',
    type=click.Choice(['right', 'bottom', 'none']),
    default='right',
    help='Legend position'
)
@click.option(
    '--no-scale',
    is_flag=True,
    default=False,
    help='Hide scale bar'
)
@click.option(
    '--width',
    type=float,
    default=12.0,
    help='Figure width in inches'
)
@click.option(
    '--chrom-height',
    type=float,
    default=0.4,
    help='Chromosome bar height in inches'
)
@click.option(
    '--font-size',
    type=float,
    default=10.0,
    help='Label font size'
)
@click.option(
    '--dpi',
    type=int,
    default=300,
    help='PNG output resolution'
)
def main(
    karyotype: Path,
    segments: Path,
    colors: Optional[Path],
    out: Path,
    sort: str,
    legend: str,
    no_scale: bool,
    width: float,
    chrom_height: float,
    font_size: float,
    dpi: int
):
    """
    HybridChromoMap - Chromosome ancestry painting visualization

    Generate publication-quality chromosome ancestry painting diagrams
    for hybrid species with variable ploidy levels.

    If no colors file is provided (-c), chromosomes will be auto-colored
    with each chromosome getting a distinct color.
    """
    try:
        # Parse input files
        click.echo(f"Loading karyotype from {karyotype}...")
        karyo = parse_karyotype(karyotype)

        click.echo(f"Loading segments from {segments}...")
        parse_segments(segments, karyo)

        # Handle colors
        if colors is not None:
            click.echo(f"Loading colors from {colors}...")
            origin_dict = parse_origins(colors)

            # Validate that segment origins match colors file
            errors = validate_data(karyo, origin_dict)
            if errors:
                click.echo("Validation errors:", err=True)
                for error in errors:
                    click.echo(f"  - {error}", err=True)
                sys.exit(1)
        else:
            click.echo("No colors file provided, auto-generating colors for origins...")
            # Generate colors based on unique origins in segments
            origin_dict = generate_origin_colors(karyo)

        # Render
        click.echo(f"Rendering to {out}...")
        renderer = ChromoMapRenderer(
            karyotype=karyo,
            origins=origin_dict,
            fig_width=width,
            chrom_height=chrom_height,
            font_size=font_size
        )
        renderer.render(
            output_path=out,
            sort_by=sort,
            legend_position=legend,
            show_scale=not no_scale,
            dpi=dpi
        )

        click.echo(f"Done! Output saved to {out}")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
