from bokeh import events
from bokeh.models import CustomJS

def point_to_px(val):
    return val * 1.33333


def color_scheme(idx):
    colors = {
        "blue": "#4e75a3",
        "lightblue": "#b3c5da",
        "green": "#7d9263",
        "lightgreen" : "#c9d3be",
        "red" : "#c00000",
        "orange" : "#dd7e0e",
        "yellow" : "#dfcf04",
        "pink" : "#b55475",
        "purple" : "#7153a1",
        "black": "#000000",
        "grey": "#555555",
    }
    return colors[idx]


def style_plot(plot):
    plot.grid.grid_line_width = point_to_px(2)
    plot.outline_line_width = point_to_px(2)
    plot.axis.axis_line_width = point_to_px(2)
    plot.axis.major_tick_line_width = point_to_px(2)
    plot.axis.minor_tick_line_width = point_to_px(1)
    plot.axis.axis_label_text_color = "black"
    plot.axis.axis_label_text_font_size = "20pt"
    plot.axis.axis_label_text_font = "calibri"
    plot.axis.axis_label_text_font_style = "normal"
    plot.axis.major_label_text_color = "black"
    plot.axis.major_label_text_font_size = "20pt"
    plot.axis.major_label_text_font = "calibri"
    plot.axis.major_tick_out = 10
    plot.axis.major_tick_in = 4
    plot.axis.minor_tick_out = 7
    plot.legend.label_text_font = "calibri"
    plot.legend.label_text_font_size = "20pt"
    plot.legend.label_text_color = "black"
    plot.legend.label_height = 0
    plot.legend.label_text_line_height = 0
    plot.legend.spacing = 0
    plot.legend.background_fill_alpha = 1
    plot.legend.border_line_alpha = 1
    plot.legend.border_line_color = "black"
    plot.legend.border_line_width = point_to_px(1)
    plot.title.text_font_size  = "20pt"
    plot.title.text_color = "black"
    plot.title.text_font = "calibri"
    plot.js_on_event(events.DoubleTap, CustomJS(args=dict(other=plot.legend[0]),
             code="other.visible = !other.visible;"
    ))

    # plot.plot_height = 365 # half high plots
