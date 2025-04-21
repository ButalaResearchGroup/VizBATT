#MIT License
# 
# Copyright (c) 2025 Danielle N. Alverson, Eric Fonseca, Kausturi Parui, Steph J. Meikle
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import dash
from dash import dcc, html, callback_context
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import wrpoly as wrp
import os
import dash_ag_grid as dag
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import io
import base64
import simpleml as sml


#os.chdir(r'C:\Users\dalverson\Dropbox (UFL)\Butala Hackathon Team\Danielles Code Folder\wrpoly\data')

merged_df = pd.read_csv('merged_new_data.csv', index_col=0)
#merged_df['Octahedral Distortion?'] = merged_df['Octahedral Distortion?'].map({'No': 0, 'Yes': 1})
unique_df = merged_df.drop_duplicates(subset='cif_name', keep='first').copy()
#unique_df['Octahedral Distortion?'] = unique_df['Octahedral Distortion?'].map({'No': 0, 'Yes': 1})

other_new_names = {'bond_angle_variance': 'Bond Angle Variance','Structure type': 'Structure Type','Average Bond Angles': 'Average Bond Length Angles','Octahedral Distortion?': 'Octahedral Distortion'}
new_column_names = {'avg_bond_length':'Average Bond Length', 'std_bond_length': 'Standard Deviation Bond Length', 'skew_bond_length': 'Skew Bond Length',
       'distortion_index': 'Distortion Index', 'quadratic_elongation': 'Quadratic Elongation', 'n_corner_pairs': 'Number of Corner Pairs',
       'n_edge_pairs': 'Number of Edge Pairs', 'bond_angle_variance': 'Bond Angle Variance', 'n_face_pairs': 'Number of Face Pairs', 'n_atoms_in_cell': 'Number of Atoms in Unit Cell', 'avg_bond_angles': 'Average Bond Angles',
       'std_bond_length_angles': 'Standard Deviation Bond Length Angles', 'skew_bond_length_angles': 'Skew Bond Length Angles', 'volumes': 'Volumes', 'central_atoms': 'Central Atoms', 'polyhedra_formula': 'Polyhedra Formula',
       'poly_types': 'Polyhedra Types', 'n_atoms_per_cell': 'Number of Atoms per Unit Cell', 'cif_name': 'Cif Name'}
merged_df.rename(columns=other_new_names, inplace=True)
unique_df.rename(columns=other_new_names, inplace=True)
merged_df.rename(columns=new_column_names, inplace=True)
unique_df.rename(columns=new_column_names, inplace=True)

merged_df['None'] = 'None'
unique_df['None'] = 'None'

unique_df.loc[:, 'for_none'] = 50
merged_df.loc[:, 'for_none'] = 50

for_none_df_u = unique_df[['for_none']]
for_none_df_m = merged_df[['for_none']]
merged_df = merged_df.drop(columns=['for_none'])
unique_df = unique_df.drop(columns=['for_none'])

unique_df.columns.metadata = {'Average Bond Length': {'description': 'Average bond length of the material', 'units': 'Angstroms'},
                              'Standard Deviation Bond Length': {'description': 'Standard deviation of bond length of the material', 'units': 'Angstroms'},
                                'Skew Bond Length': {'description': 'Skewness of bond length of the material', 'units': 'Angstroms'},
                                'Distortion Index': {'description': 'the distortion based on bond lengths from the central atom to the coordinating atom', 'units': 'Unitless'},
                                'Quadratic Elongation': {'description': 'the quadratic elongation of the octahedron', 'units': 'Unitless'},
                                'Number of Corner Pairs': {'description': 'the number of corner pairs in the octahedron', 'units': 'Unitless'},
                                'Number of Face Pairs': {'description': 'the number of face pairs in the octahedron', 'units': 'Unitless'},
                                'Number of Edge Pairs': {'description': 'the number of edge pairs in the octahedron', 'units': 'Unitless'},
                                'Number of Atoms in Unit Cell': {'description': 'the number of atoms in the unit cell', 'units': 'Unitless'},
                                'Average Bond Length Angles': {'description': 'the average bond angles of the material', 'units': 'Degrees'},
                                'Standard Deviation Bond Length Angles': {'description': 'the standard deviation of bond angles of the material', 'units': 'Degrees'},
                                'Skew Bond Length Angles': {'description': 'the skewness of bond angles of the material', 'units': 'Degrees'},
                                'Volumes': {'description': 'the volume of the unit cell', 'units': 'Angstroms^3'},
                                'Polyhedra Types': {'description': 'the polyhedral type of the material', 'units': 'Unitless'},
                                'Formula': {'description': 'the chemical formula of the material', 'units': 'Unitless'},
                                'Cif Name': {'description': 'the name of the cif file', 'units': 'Unitless'},
                                'Total Number of TM': {'description': 'the total number of transition metals in the material', 'units': 'Unitless'},
                                'Li/TM ': {'description': 'the ratio of lithium to transition metals in the material', 'units': 'Unitless'},
                                'Capacity (mAh/g)': {'description': 'the capacity of the material', 'units': 'mAh/g'},
                                '1st Charge (mAh/g)': {'description': 'the first charge of the material', 'units': 'mAh/g'},
                                '1st Discharge (mAh/g)': {'description': 'the first discharge of the material', 'units': 'mAh/g'},
                                'Octahedral Distortion': {'description': 'the octahedral distortion of the material', 'units': 'Unitless'},
                                'Central Atoms': {'description': 'the central atoms of the material', 'units': 'Unitless'},
                                'Bond Angle Variance': {'description': 'the variance of the bond angles of the material', 'units': 'Degrees'},
                                'Polyhedra Formula': {'description': 'the polyhedral formula of the material', 'units': 'Unitless'},
                                'Number of Atoms per Unit Cell': {'description': 'the number of atoms per unit cell', 'units': 'Unitless'},
                                'Mol wt (g/mol)': {'description': 'the molecular weight of the material', 'units': 'g/mol'},
                                'Central Atom(s)': {'description': 'the central atom(s) of the material', 'units': 'Unitless'},
                                '2nd Discharge (mAh/g)': {'description': 'the second discharge of the material', 'units': 'mAh/g'},
                                'x moles of Li+': {'description': 'the number of moles of lithium in the material', 'units': 'moles'},
                                'Rate of Cycling (nC)': {'description': 'the rate of cycling of the material', 'units': 'nC'},
                                'Voltage (V)': {'description': 'the voltage of the material', 'units': 'V'},
                                'Columbic Efficiency (%)': {'description': 'the columbic efficiency of the material', 'units': '%'},
                                'Reversible/Irreversible cycling': {'description': 'the reversible/irreversible cycling of the material', 'units': 'Unitless'}
}
merged_df.columns.metadata = {'Average Bond Length': {'description': 'Average bond length of the material', 'units': 'Angstroms'},
                              'Standard Deviation Bond Length': {'description': 'Standard deviation of bond length of the material', 'units': 'Angstroms'},
                                'Skew Bond Length': {'description': 'Skewness of bond length of the material', 'units': 'Angstroms'},
                                'Distortion Index': {'description': 'the distortion based on bond lengths from the central atom to the coordinating atom', 'units': 'Unitless'},
                                'Quadratic Elongation': {'description': 'the quadratic elongation of the octahedron', 'units': 'Unitless'},
                                'Number of Corner Pairs': {'description': 'the number of corner pairs in the octahedron', 'units': 'Unitless'},
                                'Number of Face Pairs': {'description': 'the number of face pairs in the octahedron', 'units': 'Unitless'},
                                'Number of Edge Pairs': {'description': 'the number of edge pairs in the octahedron', 'units': 'Unitless'},
                                'Number of Atoms in Unit Cell': {'description': 'the number of atoms in the unit cell', 'units': 'Unitless'},
                                'Average Bond Length Angles': {'description': 'the average bond angles of the material', 'units': 'Degrees'},
                                'Standard Deviation Bond Length Angles': {'description': 'the standard deviation of bond angles of the material', 'units': 'Degrees'},
                                'Skew Bond Length Angles': {'description': 'the skewness of bond angles of the material', 'units': 'Degrees'},
                                'Volumes': {'description': 'the volume of the unit cell', 'units': 'Angstroms^3'},
                                'Polyhedra Types': {'description': 'the polyhedral type of the material', 'units': 'Unitless'},
                                'Formula': {'description': 'the chemical formula of the material', 'units': 'Unitless'},
                                'Cif Name': {'description': 'the name of the cif file', 'units': 'Unitless'},
                                'Total Number of TM': {'description': 'the total number of transition metals in the material', 'units': 'Unitless'},
                                'Li/TM ': {'description': 'the ratio of lithium to transition metals in the material', 'units': 'Unitless'},
                                'Capacity (mAh/g)': {'description': 'the capacity of the material', 'units': 'mAh/g'},
                                '1st Charge (mAh/g)': {'description': 'the first charge of the material', 'units': 'mAh/g'},
                                '1st Discharge (mAh/g)': {'description': 'the first discharge of the material', 'units': 'mAh/g'},
                                'Octahedral Distortion': {'description': 'the octahedral distortion of the material', 'units': 'Unitless'},
                                'Central Atoms': {'description': 'the central atoms of the material', 'units': 'Unitless'},
                                'Bond Angle Variance': {'description': 'the variance of the bond angles of the material', 'units': 'Degrees'},
                                'Polyhedra Formula': {'description': 'the polyhedral formula of the material', 'units': 'Unitless'},
                                'Number of Atoms per Unit Cell': {'description': 'the number of atoms per unit cell', 'units': 'Unitless'},
                                'Mol wt (g/mol)': {'description': 'the molecular weight of the material', 'units': 'g/mol'},
                                'Central Atom(s)': {'description': 'the central atom(s) of the material', 'units': 'Unitless'},
                                '2nd Discharge (mAh/g)': {'description': 'the second discharge of the material', 'units': 'mAh/g'},
                                'x moles of Li+': {'description': 'the number of moles of lithium in the material', 'units': 'moles'},
                                'Rate of Cycling (nC)': {'description': 'the rate of cycling of the material', 'units': 'nC'},
                                'Voltage (V)': {'description': 'the voltage of the material', 'units': 'V'},
                                'Columbic Efficiency (%)': {'description': 'the columbic efficiency of the material', 'units': '%'},
                                'Reversible/Irreversible cycling': {'description': 'the reversible/irreversible cycling of the material', 'units': 'Unitless'}
}
                                
image_path ='NSF_Official_logo_High_Res_1200ppi.png'                               
# Using base64 encoding and decoding
def b64_image(image_filename):
    with open(image_filename, 'rb') as f:
        image = f.read()
    return 'data:image/png;base64,' + base64.b64encode(image).decode('utf-8')
    

# Assume 'unique_df' is your dataframe with data
# Create Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])


# Options for dropdowns
options = [{'label': col, 'value': col} for col in merged_df.columns]

# Define the layout of the app
app.layout = dbc.Container([
    
    html.Div([
        html.Div([
            html.A(html.H3("Butala Research Group", style={'color': 'white', 'margin-left': '40px'}),href="https://butala.mse.ufl.edu/", target="_blank", style={'color': 'white'}),
        ], style={'display': 'flex', 'flex-direction': 'column','margin-top': '20px'}),
        html.Div([
            dbc.Label("Polyhedra Analysis Software", style={'color': 'white', 'margin-left': '30px', 'margin-right': '20px'}),
            html.A(html.H3("WRPoly", style={'color': 'white', 'margin-left': '30px', 'margin-right': '20px'}), href="https://github.com/ericfonseca95/wrpoly", target="_blank", style={'color': 'white'}),
        ], style={'margin-top': '10px'}),
    ], style={'background-color': '#1d7a72', 'width': '100%', 'display': 'flex', 'justify-content': 'space-between'}),
    html.Br(),  # Add a line break

    html.H1("Battery Materials Dashboard", style={'text-align': 'center'}),
    
    # dbc.Row([
    #     dbc.Col([
    #         dcc.Tabs(id='tabs', value='tab-1', children=[
    #             dcc.Tab(label='Scatter Plot', value='tab-1'),
    #             dcc.Tab(label='Linear Regression', value='tab-2')
    #         ])
    #     ], width={'size': 6, 'offset': 3})  # Adjust width and offset as needed
    # ]),
    
    html.Div(id='tabs-content'),
    
    dbc.Row([
        dbc.Col([
            dbc.Label("Select Displayed Material Type"),
            dcc.Dropdown(
                id='material-type-dropdown',
                options=[ {'label': 'Wadsley-Roth Materials', 'value': 'Wadsley-Roth'},{'label': 'Wadsley-Roth Adjacent Materials', 'value': 'Wadsley-Roth Adjacent'}, {'label': 'All Materials', 'value': 'all'}],
                value = 'all',
                style = {'width': '100%', 'background-color': '#a5cac7', 'color': 'black'},
            ),
            html.Br(),  # Add a line break
            
            dbc.Label("Select X Parameter:"),
            dcc.Dropdown(
                id='x-dropdown',
                options=options,
                value='Distortion Index',  # Default selected value
                style={'width':'100%', 'background-color': '#a5cac7', 'color': 'black'},  # Set width of the dropdown
            ),
            html.Br(),  # Add a line break
            
            dbc.Label("Select Y Parameter:"),
            dcc.Dropdown(
                id='y-dropdown',
                options=options,
                value='Capacity (mAh/g)',  # Default selected value
                style={'width': '100%', 'background-color': '#a5cac7', 'color': 'black'},
            ),
            html.Br(),  # Add a line break
            
            dbc.Label("Select Size Parameter:"),
            dcc.Dropdown(
                id='size-dropdown',
                options=options,
                value='Volumes',  # Default selected value
                style={'width': '100%','background-color': '#a5cac7', 'color': 'black'},
            ),
            html.Br(),  # Add a line break
            
            dbc.Label("Select Color Parameter:"),
            dcc.Dropdown(
                id='color-dropdown',
                options=options,
                value='Formula',  # Default selected value
                style={'width': '100%', 'background-color': '#a5cac7', 'color': 'black'},
            ),
            html.Br(),  # Add a line break
            
            dbc.Label("Select Shape Parameter:"),
            dcc.Dropdown(
                id='shape-dropdown',
                options=options,
                value='Polyhedra Types',  # Default selected value
                style={'width': '100%', 'background-color': '#a5cac7', 'color': 'black'},
            ),
            html.Br(),  # Add a line break
            
            dbc.Label("Show Legend:"),
            dbc.Checklist(
                id='toggle',
                options=[
                    {'label': 'Yes', 'value': True},
                ],
                value=[False], # Default to 'No'
                switch = True, 
            ),
            html.Br(),  # Add a line break
            
            dbc.Label("Search Formula:"),
            dbc.Input(id='search-bar', 
                      type='text', 
                      placeholder='Enter formula'),
            html.Br(),
            
            html.Div(id='upload-status'), #Display upload status
        
            dbc.Label("Upload CIF File:"),
            dcc.Upload(id='upload-data',
                children=html.Div([
                'Drag and Drop or ',
            html.A('Select a CIF File')
                ]),
            style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
            multiple=False
            ),
            
        html.Br(),  # Add a line break
            
            # Button here
        dbc.Label("Click ONCE after importing cif file", style={'color': '#fc6c85', "font-size": "15px", 'margin-left': '20px', 'margin-right': '20px'}),
        dbc.Button("I believe in fairies", color="pink", style={"font-size": "20px", "color": "pink", "border-radius": "20px", "box-shadow": "0px 0px 15px pink", "width": "100%", "margin-bottom": "10px"}),
            
            
        ], width = 2),  # Set to True if you want to allow multiple file uploads
        
        # Scatter Plot
        dbc.Col([
            dcc.Graph(id='scatter-plot')
        ], width=10),
        
    ], className="mb-4", style = {'background-color': '#e9eaec', 'margin': '20px','padding': '10 px', 'border-radius': '10px', 'box-shadow': '0px 0px 15px #e9eaec'}),
     
    html.Br(),
        # Add bottom border
        

    html.Div(style={'border-top': '2px solid black', 'margin-top': '20px', 'padding-top': '20px', 'text-align': 'center'}, children=[
    #calling the above function
        html.Img(src=b64_image(image_path),style={'text-align':'center','height':'10%', 'width':'10%'}),
        html.P("NSF DMR-2334240", style={'text-align': 'center'}),
        html.P("Authors: Danielle N. Alverson, Kausturi Parui, Eric Fonseca Steph J. Meikle, Megan M. Butala", style={'text-align': 'center'}),
        html.P("This project is licensed under the MIT License.", style={'text-align': 'center'}),
    ]),
    # new buttons for ml
    ],
    
 fluid=True)


@app.callback(
    [Output('x-dropdown', 'options'),
     Output('y-dropdown', 'options'),
     Output('size-dropdown', 'options'),
     Output('color-dropdown', 'options'),
     Output('shape-dropdown', 'options')],
     [Input('x-dropdown', 'value'),
     Input('y-dropdown', 'value'),
     Input('size-dropdown', 'value'),
     Input('color-dropdown', 'value'),
     Input('shape-dropdown', 'value')],
)
def perform_action_and_update_options(selected_x, selected_y, selected_size, selected_color, selected_shape):
    # Generate options for dropdowns
    options_with_tooltips = []
    for col in unique_df.columns:
        label = col
        value = col
        metadata = unique_df.columns.get_level_values(0).metadata.get(col, {})
        title = metadata.get('description', f'This selected parameter {col}')
        option = {'label': label, 'value': value, 'title': title}
        options_with_tooltips.append(option)

    return (options_with_tooltips,
            options_with_tooltips,
            options_with_tooltips,
            options_with_tooltips,
            options_with_tooltips)


# Callback to update scatter plot based on dropdown selection
@app.callback(
    [Output('scatter-plot', 'figure'),
     Output('upload-data', 'contents'),
     Output('upload-status', 'children')],
    [Input('material-type-dropdown', 'value'),
     Input('x-dropdown', 'value'),
     Input('y-dropdown', 'value'),
     Input('size-dropdown', 'value'),
     Input('color-dropdown', 'value'),
     Input('shape-dropdown', 'value'),
     Input('toggle', 'value'),
     Input('search-bar', 'value')],
    [State('upload-data', 'contents'),
     State('upload-data', 'filename')]
)

def update_scatter_plot(wr_mat, selected_x, selected_y, selected_size, color, symbol, toggle_value, search_term, file_contents, filename):
    print(wr_mat)

        
    # new_column_names = {'avg_bond_length':'Average Bond Length', 'std_bond_length': 'Standard Deviation Bond Length', 'skew_bond_length': 'Skew Bond Length',
    #    'distortion_index': 'Distortion Index', 'quadratic_elongation': 'Quadratic Elongation', 'n_corner_pairs': 'Number of Corner Pairs',
    #    'n_edge_pairs': 'Number of Edge Pairs', 'bond_angle_variance': 'Bond Angle Variance', 'n_face_pairs': 'Number of Face Pairs', 'n_atoms_in_cell': 'Number of Atoms in Unit Cell', 'avg_bond_angles': 'Average Bond Angles',
    #    'std_bond_length_angles': 'Standard Deviation Bond Length Angles', 'skew_bond_length_angles': 'Skew Bond Length Angles', 'volumes': 'Volumes', 'central_atoms': 'Central Atoms', 'polyhedra_formula': 'Polyhedra Formula',
    #    'poly_types': 'Polyhedra Types', 'n_atoms_per_cell': 'Number of Atoms per Unit Cell', 'cif_name': 'Cif Name'}
    
    if wr_mat  == 'all':
        filtered_df = merged_df
    if wr_mat == 'Wadsley-Roth':
        filtered_df = merged_df[(merged_df['Structure'] == 'Wadsley-Roth')]
        
    if wr_mat == 'Wadsley-Roth Adjacent':
        filtered_df = merged_df[merged_df['Structure'] == 'Wadsley-Roth Adjacent']

    # Reset index to ensure unique indices
    filtered_df.reset_index(drop=True, inplace=True)

    # Print duplicated rows based on index
    #print(filtered_df[filtered_df.index.duplicated()])

    # Ensure unique columns by dropping duplicates
    filtered_df = filtered_df.loc[:, ~filtered_df.columns.duplicated()]

      # Check if selected_size is None or 'None'
    if selected_size is None or selected_size == 'None':
        selected_size = for_none_df_u['for_none']  # Or any other default size value you want to use
    
    show_legend = True in toggle_value if toggle_value else False
    
    # Check if search_term is None or empty
    if search_term is None or search_term == '':
        filtered_df_2 = filtered_df  # Display all data
        print(search_term)
    else:
        # Filter data based on search term
        filtered_df_2 = merged_df[(merged_df['Cif Name'].str.contains(search_term, case=False, na=False)) |
                        (merged_df['Formula'].str.contains(search_term, case=False, na=False))]
        print(search_term)
        
    # Check if file_contents is None
    if file_contents is None:
        print(file_contents)
    else:
        content_type, content_string = file_contents.split(',')
        decoded = base64.b64decode(content_string)
        file_extension = filename.split('.')[-1]
        cif_names = filename.split('.')[0]
        print(file_extension)
        print(cif_names)
        
        if file_extension == 'cif':
          # Read the CIF file
            temp_filepath = r'C:\Users\dalverson\Dropbox (UFL)\Butala Hackathon Team\Danielles Code Folder\wrpoly\temp\temp.cif'
            if os.path.exists(temp_filepath):
                os.remove(temp_filepath)
                file_contents = None
                print('File already exists')
            else:
                
                with open(temp_filepath, 'wb') as temp_file:
                    temp_file.write(decoded)
                
                 # Read the CIF file and get the structure
                file = io.StringIO(decoded.decode('utf-8'))
                
    #             # Add appropriate import statements for these functions
    # get_cif_name(file)
                #cif_names = get_cif_name(file)
                #print(cif_names)
                structure = wrp.get_structure(temp_filepath)
                
                uploaded_df = wrp.get_average_df(structure, cif_names)
                print(uploaded_df)
                
    #             # Add appropriate definition for new_column_names
                uploaded_df.rename(columns=new_column_names, inplace=True)
                
    #             # Merge the uploaded data with the filtered data
                if filtered_df_2 is None:
                    filtered_df_2 = pd.concat([merged_df, uploaded_df], ignore_index=True)
                else:
                    filtered_df_2 = pd.concat([filtered_df, uploaded_df], ignore_index=True)
                os.remove(temp_filepath)
                file_contents = None
                
               # if selected_size == 'None':
                    #selected_size = 100
         
    fig = px.scatter(filtered_df_2, x=selected_x, y=selected_y, hover_data=['Cif Name'], size=selected_size, color=color, symbol=symbol,
                     color_continuous_scale='Jet')
    
    #fig = px.scatter(merged_df, x=selected_x, y=selected_y, hover_data=['cif_name'], size=size, color=color, symbol=symbol,
                     #color_continuous_scale='Jet')
    
    # Customizing the font
    fig.update_layout(font_family="Helvetica", font_color="black", title_font_family="Helvetica", title_font_color="black")
    
    # Customizing theme
    fig.layout.template = 'presentation'
    
    # Customizing the size of the plot
    fig.update_layout(width=1500, height=1000)
    
    # Customize title
    if selected_size is not None:
        fig.update_layout(title=(selected_x + ' vs. '+selected_y), title_x=0.5)
        

    # Customize legend
    fig.update_layout(legend_title_text=None)
    if show_legend is not None:
        fig.update_layout(showlegend=show_legend)

    # Customizing margin
    fig.update_layout(margin=dict(l=70, r=70, t=70, b=70))
    
    
    
    return fig, file_contents, f'File "{filename}" uploaded and processed.'

def get_cif_name(filename):
    if isinstance(filename, io.StringIO):
        # Handle StringIO differently
        return filename.split('.cif')[0] # You may need to replace this with your logic for StringIO
    else:
        # Assuming filename is a string representing a file path
        return filename.split('/')[-1].split('.cif')[0]

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, port=8057)