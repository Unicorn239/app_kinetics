# Unicorn
# Created for kinetics model building and prediction

import os
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import pandas as pd
import numpy as np
import re
from scipy.optimize import curve_fit
import plotly.graph_objs as go


app = dash.Dash(__name__)
server = app.server
port = int(os.environ.get("PORT", 5000))

app.layout = html.Div(
             children = [
             html.H3('Enter reaction time here...'),
             dcc.Textarea(id = 'time',
                          placeholder = 'Enter reaction time',
                          value = '',
                          style = {
                                   'width' : '45%',
                                   'height' : 30
                                      }
                             ),
             html.H3('Enter conversion here...'),
             dcc.Textarea(id = 'conv',
                       placeholder = 'Enter conversion',
                       value = '',
                       style = {
                                 'width' : '45%',
                                  'height' : 30
                                      }                       
                             ),
             html.H3('Parameters after fitting are shown as below...'),
             html.Button('calculate', id = 'button'),
             html.Div(id = 'results', 
                      style = {'fontSize' : '16x'}),
             html.H3('Conversion curve is shown as below...'),      
             dcc.Graph(id = 'curve')
                        ])


@app.callback(
             Output('results', 'children'),
             Input('time', 'value'),
             Input('conv', 'value'),
             Input('button', 'n_clicks')
              )
def update_params(input_time, input_conv, n_clicks):
    mo = re.compile(r'\d*\.?\d+')
    lst_time = mo.findall(input_time)
    lst_conv = mo.findall(input_conv)
    time_raw = pd.Series(lst_time).astype(float)
    conv_raw = pd.Series(lst_conv).astype(float)
    df = pd.DataFrame(
                      {'time': time_raw,
                       'conv': conv_raw}
                      )
    def fst_order(t, Cmax, k, rsd):
        return Cmax*(1-np.exp(-k*t)) + rsd
    params1, cov1 = curve_fit(fst_order, df.time, df.conv, p0 = [100, 0.1,\
                                                                 0])
    SS_tot1 = ((df.conv - df.conv.mean())**2).sum() 
    SS_rsd1 = ((fst_order(df.time, params1[0], params1[1], params1[2])-df.\
               conv)**2).sum()
    r_square1 = 1 - SS_rsd1/SS_tot1 
    
    def sec_order(t, Cmax, k, rsd):
        return Cmax**2*k*t/(1+Cmax*k*t) + rsd
    
    params2, cov2 = curve_fit(sec_order, df.time, df.conv, p0 = [100, 0.1,\
                                                                 0])
    SS_tot2 = ((df.conv - df.conv.mean())**2).sum() 
    SS_rsd2 = ((sec_order(df.time, params2[0], params2[1], params2[2])-df.\
               conv)**2).sum()
    r_square2 = 1 - SS_rsd2/SS_tot2
    
    if n_clicks > 0:
        if r_square1 > r_square2 and r_square1 > 0.95:
            return ['first-order kinetics', '\t', 'Cmax = ', params1[0].round(3),\
                    '\t', 'k = ', params1[1].round(3), '\t', 'rsd = ', \
                    params1[2].round(3), '\t', 'R^2 = ', r_square1.round(3)]
        elif r_square2 > r_square1 and r_square2 > 0.95:
            return ['second-order kinetics', '\t', 'Cmax = ', params2[0].round(3),\
                    '\t', 'k = ', params2[1].round(3), '\t', 'rsd = ', \
                        params2[2].round(3), '\t', 'R^2 = ', r_square2.round(3)]
        elif r_square1 <= 0.95 and r_square2 <= 0.95:
            return 'no good kinetic curve fitting, contact chemical engineer'
        

@app.callback(
              Output('curve', 'figure'),
              [Input('time', 'value'),
               Input('conv', 'value')]
              )
def update_figure(input_time, input_conv):
    mo = re.compile(r'\d*\.?\d+')
    lst_time = mo.findall(input_time)
    lst_conv = mo.findall(input_conv)
    time_raw = pd.Series(lst_time).astype(float)
    conv_raw = pd.Series(lst_conv).astype(float)
    df = pd.DataFrame(
                      {'time': time_raw,
                       'conv': conv_raw}
                      )
    trace1 = go.Scatter(x = df.time, y = df.conv, mode = 'markers',
                        marker = {'size' : 5})
    def fst_order(t, Cmax, k, rsd):
        return Cmax*(np.exp(-k*t)) + rsd
    params1, cov1 = curve_fit(fst_order, df.time, df.conv, p0 = [100, 0.1,\
                                                                 0])
    SS_tot1 = ((df.conv - df.conv.mean())**2).sum() 
    SS_rsd1 = ((fst_order(df.time, params1[0], params1[1], params1[2])-df.\
               conv)**2).sum()
    r_square1 = 1 - SS_rsd1/SS_tot1 
    
    def sec_order(t, Cmax, k, rsd):
        return Cmax**2*k*t/(1+Cmax*k*t) + rsd
    
    params2, cov2 = curve_fit(sec_order, df.time, df.conv, p0 = [100, 0.1,\
                                                                 0])
    SS_tot2 = ((df.conv - df.conv.mean())**2).sum() 
    SS_rsd2 = ((sec_order(df.time, params2[0], params2[1], params2[2])-df.\
               conv)**2).sum()
    r_square2 = 1 - SS_rsd2/SS_tot2
    
    interval = np.linspace(df.time.min(), df.time.max())
    
    if r_square1 > r_square2 and r_square1 > 0.95:
        conv_fit = fst_order(interval, params1[0], params1[1], params1[2])
        df_fit = pd.DataFrame({'time': interval, 'conv': conv_fit})       
        trace2 = go.Scatter(x = df_fit.time, y = df_fit.conv, mode = 'lines',
                        line = {'width' : 1})          
        fig = go.Figure([trace1, trace2])
        return fig
    elif r_square2 > r_square1 and r_square2 > 0.95:
        conv_fit = sec_order(interval, params2[0], params2[1], params2[2])
        df_fit = pd.DataFrame({'time': interval, 'conv': conv_fit})       
        trace2 = go.Scatter(x = df_fit.time, y = df_fit.conv, mode = 'lines',
                        line = {'width' : 1})          
        fig = go.Figure([trace1, trace2])
        return fig
    else:
        fig = go.Figure(trace1)
        return fig
            

if __name__ == '__main__':
    kinetics.run_server(debug = False, 
                   host="0.0.0.0",
                   port=port)
    
