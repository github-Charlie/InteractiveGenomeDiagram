from __future__ import print_function

from bokeh.document import Document
from bokeh.embed import file_html
from bokeh.resources import INLINE
from bokeh.browserlib import view
from bokeh.models.glyphs import Arc, Text, AnnularWedge, Line
from bokeh.models import ColumnDataSource, Range1d, Plot, PanTool, WheelZoomTool, PreviewSaveTool, ResetTool, \
    HelpTool, HoverTool

from Bio import SeqIO
from math import pi, sin, cos


def polar_to_cartesian(r, alpha):
    return r*cos(alpha), r*sin(alpha)

xdr = Range1d(start=-1, end=1)
ydr = Range1d(start=-1, end=1)

plot = Plot(title="test diagram", x_range=xdr, y_range=ydr, plot_width=800, plot_height=800)

record = SeqIO.read('AF488832.gbk', 'gb')

track_mid_radius = 0.7
track_scale = 0.1
track_inner_radius = track_mid_radius - track_scale
track_outer_radius = track_mid_radius + track_scale

start_angles = []
end_angles = []
wedge_inner_radii = []
wedge_outer_radii = []

type_list = []
start_locations = []
end_locations = []
strand_values = []
locus_tags = []

for feature in record.features:
    # Ignore source features and features spanning the whole record.
    if feature.type == 'source':
        continue
    if feature.location.start == 0 and feature.location.end == len(record):
        continue

    start_angle = feature.location.start/len(record) * 2*pi
    end_angle   = feature.location.end  /len(record) * 2*pi
    # Bokeh measures glyph angles anti-clockwise from horizontal (3 O'clock) for some reason.
    # I want to measure clockwise from vertical, so have to apply a couple of transformations:
    # take the negative of the angle and add 90 degrees
    start_angle = pi/2 - start_angle
    end_angle   = pi/2 - end_angle

    if feature.location.strand == -1:           # Put these on the inside
        wedge_inner_radius = track_inner_radius
        wedge_outer_radius = track_mid_radius
        text_align = "right"
        text_radius = wedge_inner_radius - 0.01
    else:                                       # And these on the outside
        wedge_inner_radius = track_mid_radius
        wedge_outer_radius = track_outer_radius
        text_align = "left"
        text_radius = wedge_outer_radius + 0.01

    if len(locus_tags) >= 1:
        if locus_tags[-1] == feature.qualifiers['locus_tag'][0]:  # If there is already a feature with the tag name
            type_list[-1] = feature.type                          # Override it
            continue

    # Due to a tooltip bug, it was not possible to generate the AnnularWedge objects one by one as I have done
    # with the text.  Hence the many lists to pass to the ColumnDataSource.
    wedge_inner_radii.append(wedge_inner_radius)
    wedge_outer_radii.append(wedge_outer_radius)
    start_angles.append(start_angle)
    end_angles.append(end_angle)

    type_list.append(feature.type)
    start_locations.append(feature.location.start)
    end_locations.append(feature.location.end)
    strand_values.append(feature.strand)
    locus_tags.append(feature.qualifiers['locus_tag'][0])

    #Add text labels
    x,y = polar_to_cartesian(r=text_radius, alpha=end_angle)
    text = Text(x=x,
                y=y,
                text=[feature.qualifiers['locus_tag'][0]],
                angle=end_angle,
                text_font_size='10px',
                text_align=text_align)
    plot.add_glyph(text)

source = ColumnDataSource(dict(type=type_list,
                               start_location=start_locations,
                               end_location=end_locations,
                               strand=strand_values,
                               locus_tag=locus_tags,
                               start_angle=end_angles,
                               end_angle=start_angles,
                               wedge_inner_radius=wedge_inner_radii,
                               wedge_outer_radius=wedge_outer_radii))

annular_wedge = AnnularWedge(x=0, y=0,
                             inner_radius='wedge_inner_radius',  # For the inner_radius, use the values in the column
                                                                 # called 'wedge_inner_radius'
                             outer_radius='wedge_outer_radius',  # Likewise for the outer radius...
                             line_color='white',
                             fill_color='green',
                             start_angle='start_angle',          # ...and start and end angles
                             end_angle='end_angle',
                             direction='anticlock')

features = plot.add_glyph(source, annular_wedge)  # Add all the wedges to the plot

# Add track
arc = Arc(x=0, y=0,
          radius=track_mid_radius,
          start_angle=0,
          end_angle=2*pi,
          direction='clock',
          line_width=2,
          line_color="black")
plot.add_glyph(arc)

# Add start tick
source = ColumnDataSource(dict(x=[0, 0], y=[track_inner_radius, track_outer_radius]))
line = Line(x='x', y='y', line_color="black", line_width=2)
plot.add_glyph(source, line)

# Add tooltip
features_tooltip = [
    ("Start base", "@start_location"),
    ("End base", "@end_location"),
    ("Type", "@type"),
    ("Locus tag", "@locus_tag"),
    ("Strand", "@strand")
]

# Set renderers=[features] so that the tooltip appears only for Annular Wedge glpyhs
features_hovertool = HoverTool(renderers=[features], tooltips=features_tooltip, point_policy='follow_mouse')

#Add tools to the plot; create a document; add the plot to it
plot.add_tools(WheelZoomTool(), PreviewSaveTool(), ResetTool(), HelpTool(), PanTool(), features_hovertool)
doc = Document()
doc.add(plot)


if __name__ == "__main__":
    filename = "test.html"
    with open(filename, "w") as f:
        f.write(file_html(doc, INLINE, "test"))
    print("Wrote %s" % filename)
    view(filename)

