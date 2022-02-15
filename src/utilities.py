from dash import dcc, html
import dash_bootstrap_components as dbc


def create_parameter_sliders(parameters):
    """Generate a list of rows containing a number display and slider for each
    parameter

    Parameters
    ----------
    parameters : dict
        Nested dictionary with parameter names as keys and a dict for value.
        The value-dict has keys determining the initial value, min, and max.
        The keys of the value-dict have to be: 'init', 'min', and 'max'.
        For example:
        parameters = {'param_a' : {'init': 1.0, 'min': 0.0, 'max': 2.0}}

    Returns
    -------
    list
        A list of dbc.Row() elements. Each row contains a display that shows
        the current parameter value and a slider. The ID of the display is
        '<param_name>-display'. Similarly, the ID of the slider is
        '<param_name>-slider'.
    """

    all_sliders = []
    for name in parameters.keys():
        init_val = parameters[name]['init']
        min_val = parameters[name]['min']
        max_val = parameters[name]['max']
        step = (max_val - min_val) / 1000
        display = html.H5(f"{name}: {init_val:.3f}", id=f"{name}-display",
                          style={"padding-left": "10%"})
        slider = dcc.Slider(id=f'{name}-slider', min=min_val, max=max_val,
                            step=step, value=init_val)
        all_sliders.append(dbc.Row([display, slider], style={"width": "100%"}))
    return all_sliders


class Parameters:
    """Object that stores a parameter set of an ODE"""

    def __init__(self, p_dict):
        """
        Parameters
        ----------
        p_dict : dict
            Dictionary with parameter names as keys and their values as
            dict-values. For example: p_dict = {'p_a': 1.0, 'p_b': 2.0}.
        """

        for p_name, p_val in p_dict.items():
            setattr(self, p_name, p_val)
