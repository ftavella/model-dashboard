import sys
import dash
import importlib
from dash import dcc, html
import plotly.express as px
from scipy.integrate import solve_ivp
import dash_bootstrap_components as dbc
from dash.dependencies import Output, Input
from utilities import Parameters, create_parameter_sliders

if len(sys.argv) != 2:
    raise ValueError("Please provide input file")
elif len(sys.argv) == 2:
    file_location = '/'.join(sys.argv[1].split('/')[:-1])
    sys.path.append(file_location)
    file_name = sys.argv[1].split('/')[-1]
    module_name = file_name.split('.')[0]
    try:
        model = importlib.import_module(module_name)
    except ImportError:
        raise ImportError("File not found")


colors = px.colors.qualitative.Dark24
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "Mathematical Model Interactive Dashboard"
server = app.server

app.layout = dbc.Container(
    [
        dbc.Row([
            dbc.Col([
                html.H1(app.title, style={"text-align": "center"}),
            ], className="py-2"),
        ], className="g-0", style={"flex": "0 1 auto"}),
        dbc.Row([
            dbc.Col([
                dbc.Row([
                    dbc.Col([html.H4("Model Parameters")],
                            style={"text-align": "center"},
                            className="pt-1"),
                    dbc.Col([dbc.Button('Reset', id='reset_params',
                            class_name='primary my-1', n_clicks=0)], width=3),
                ], style={"width": "100%"}, className="pb-3"),
                *create_parameter_sliders(model.parameters),
                dbc.Row([
                    dbc.Col([html.H4("Simulation Parameters")],
                            style={"text-align": "center"})
                ], style={"width": "100%"}, className="py-3"),
                dbc.Row([
                    dbc.Col([
                        dbc.Row(html.H5("Simulation Time",
                                style={"text-align": "center"}),
                                className="pt-2"),
                    ]),
                    dbc.Col([
                        dbc.Input(id="sim_time-input", type="number",
                                  value=1000, min=1e-4, max=1e4,
                                  debounce=False)
                    ]),
                ], style={"width": "100%"}, className="px-3", align="center"),

                dbc.Row([
                    html.Div([html.H4("Plot Parameters")],
                             style={"text-align": "center"})
                ], style={"width": "100%"}, className="py-3"),
                dbc.Row([
                    dbc.Col([
                        dbc.Row(html.H5("Plot height (px)",
                                style={"text-align": "center"}),
                                className="pt-2"),
                    ]),
                    dbc.Col([
                        dbc.Input(id="plot_height-input", type="number",
                                  value=200, min=10, max=800, debounce=False)
                    ]),
                ], style={"width": "100%"}, className="px-3", align="center"),
            ], width=3, style={"height": "100%", "overflow": "auto"}),

            dbc.Col([
                dbc.Row(
                    dcc.Graph(
                        id=f"{var}"), id=f"{var}-row", style={"width": "100%"})
                    for var in model.variables.keys()],
                    width=9, style={"height": "100%", "overflow": "auto"}),
        ], className="g-0", style={"flex": "1 1 auto", "overflow": "auto"}),
    ], fluid=True, className="g-0",
    style={"height": "100vh", "display": "flex", "flex-direction": "column"})


@app.callback(
    [Output(f"{p}-display", "children") for p in model.parameters.keys()],
    [Input(f"{p}-slider", "drag_value") for p in model.parameters.keys()],
    prevent_initial_call=True)
def update_labels(*params):
    return [f"{p}: {params[i]:.3f}"
            for i, p in enumerate(model.parameters.keys())]


@app.callback(
    [Output(f"{p}-slider", "value") for p in model.parameters.keys()],
    Input('reset_params', 'n_clicks'),
    prevent_initial_call=True)
def reset_params(n_clicks):
    return [model.parameters[p]['init'] for p in model.parameters.keys()]


@app.callback(
    [Output(f"{v}", "figure") for v in model.variables.keys()],
    Input('sim_time-input', "value"), Input("plot_height-input", "value"),
    [Input(f"{p}-slider", "value") for p in model.parameters.keys()],
    prevent_initial_call=False)
def simulate(t_stop, height, *params):
    if t_stop is None:
        t_stop = 1000
    elif height is None:
        height = 200

    p_dict = {p: float(val) for p, val in
              list(zip(model.parameters.keys(), params))}
    p = Parameters(p_dict)
    init_cond = [val for _, val in model.variables.items()]
    sol = solve_ivp(model.equations, [0, t_stop], init_cond, args=(p,),
                    method='LSODA')

    figures = []

    for idx, var in enumerate(model.variables.keys()):
        fig = px.scatter(x=sol.t, y=sol.y[idx, :],
                         labels={'x': 'Time', 'y': f"{var}"},
                         height=height)
        fig['data'][0]['line']['color'] = colors[idx]
        fig.update_layout(font=dict(size=15), margin=dict(t=0.1, b=0, r=25),
                          hovermode=False)
        figures.append(fig)

    return figures


if __name__ == "__main__":
    app.run_server(debug=True)
