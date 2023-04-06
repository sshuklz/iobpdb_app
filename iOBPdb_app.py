import dash
from dash import dash_table
from dash import dcc
from dash import html
from dash import MATCH, ALL
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
import dash_mantine_components as dmc
from dash.dependencies import (Input, Output, State)
import dash_bio
import dash_bio as dashbio
from dash_bio.utils import pdb_parser as parser, mol3dviewer_styles_creator as sparser

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from shutil import copy2
import pandas as pd
import os
from molmass import Formula
import urllib.request

import warnings
warnings.filterwarnings("ignore", category=UserWarning)



app=dash.Dash(__name__,suppress_callback_exceptions=True)
server=app.server

dark_theme = dmc.MantineProvider(withGlobalStyles=True, theme={"colorScheme": "dark"})



cob_df = pd.read_csv('Compound_OBP_binding.csv', dtype = str) 
cob_df['id'] = cob_df['CAS-number']
ci_df = pd.read_csv('Compound_info.csv', dtype = str) 
oi_df = pd.read_csv('OBP_info.csv', dtype = str) 

compact_comp_df = ci_df[ci_df.columns[[0,1,3,2]]]
compact_comp_df['id'] = compact_comp_df['CAS-number']
compact_comp_df['Compound name'] = compact_comp_df['Compound name'].str.split(' /').apply(lambda x: x[0])
compact_comp_df['Compound name'] = compact_comp_df['Compound name'].str.capitalize()
compact_comp_df['OBPs'] = compact_comp_df['OBPs'].astype(float)


compact_OBP_df = oi_df[oi_df.columns[[0,1,9,10,5]]]
compact_OBP_df['id'] = compact_OBP_df['Binding Protein Name']
compact_OBP_df = compact_OBP_df.sort_values('Binding Protein Name')



def unique_list(df,df_col_name):

    options = []
    
    name_list = sorted(list(set(list(df[df_col_name]))))
    
    for option in name_list:
        
        if option != '-':
        
            options += [{
                'label': option,
                'value': option
            }]
        
    return options



fg_list = ['acetylenic carbon','aldehyde','amide','amino acid','azo nitrogen',
'azole','diazene','bromine','carbamate','carbamic ester','carboxylate ion',
'carbonyl group','carbonyl with carbon','carbonyl with nitrogen',
'carbonyl with oxygen','carboxylic acid','chlorine','di sulfide','enol',
'ester','ether','fluorine','hydrazine','hydroxyl','hydroxyl in alcohol',
'hydroxyl in carboxylic acid','imine','ketone','phenol','primary alcohol',
'primary amine','mono sulfide','nitrile','nitro','secondary amine','sulfide',
'sulfonamide','sulfoxide','sulfuric acid ester','thioamide','vinylic carbon']

hp_list = ['Faeces','Urine','Breath','Skin','Milk','Blood','Saliva']

ds_list = ['SARS-CoV-2','Ubiquitous viral VOCs','rhinovirus',
           'Influenza virus H9N2, H6N2, H1N1 (H)','Respiratory syncytial virus',
           'Ubiquitous bacterial VOCs','Staphylococus Aureus',
           'Streptococus pneumoniae','Clostridium difficile',
           'Mycobaterium Tuberculosis','Haemophilus Influenzae',
           'Escherichia coli','Klebisiella pneumoniae',
           'Pseudomonas aeruginosa','COPD','Lung cancer','Asthma']

cs_list = ['3','4','5','6','7','8','9','10','11','12','19']



DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PDB_iodbp')

pdb_options = []

PDB_list = ['1KX9','3L4L','3QME','3L4A','3D73','3D74','3FE6','3FE8','3FE9',
            '3CZ0','3BFA','3BFB','2H8V','3CAB','3CDN','3BJH','3CZ1','3BFH',
            '3D75','3D76','3D77','3D78','4Z39','6QQ4','3N7H','5EL2','4FQT',
            '3K1E','1N8U','1N8V','4IJ7','5DIC','3B6X','6HHE','3OGN','6P2E',
            '6OPB','6OGH','6OMW','3R72','2WC5','3S0D','3S0E','3RZS','3S0B',
            '3B88','3B87','3B86','2QDI','2GTE','1OOF','6OTL','6OII','4PT1',
            '1OOG','3Q8I','3VB1','4F7F','3V2L','6VQ5','1OOH','4INX','4INW',
            '2P71','2P70','1DQE','1ORG','2WCH','2WCL','2WCM','2WCJ','2WC6',
            '3R1O','3R1P','1P28','1OW4','3CYZ']

for PDB in PDB_list:
    
    pdb_options += [{
        'label': PDB,
        'value': os.path.join(DATAPATH, PDB +'.pdb')
    }]

Navbar = dmc.Navbar(
    
    p=10,
    fixed = True,
    width={"base": 200},
    children=[
        
        dmc.Title("iOBPdb", order=1),
        dmc.Divider(color="red",style = {'paddingTop': 15}),
        
        dmc.Group(
            direction="column",
            children=[
                
                dmc.Anchor('Home', href='/Home'),
                dmc.Anchor('Search', href='/Search'),
                dmc.Anchor('Browse Compounds', href='/BrowseCompounds'),
                dmc.Anchor('Browse OBPs', href='/BrowseOBPs'),
                dmc.Anchor('Data Entry Form', href='/DataEntry'),
                dmc.Anchor('Download Data Files', href='/Download')
            
            ]
        )
    ],
)

app.layout = html.Div([
    dark_theme,
    
    dmc.MantineProvider(
    theme={"colorScheme": "dark"},
    children=[
        dcc.Store(id='row_selected_sel_OBP', storage_type='memory'), 
        dcc.Store(id='row_selected_sel_comp', storage_type='memory'), 
        dcc.Store(id='row_selected_OBP', storage_type='memory'), 
        dcc.Store(id='row_selected_comp', storage_type='memory'), 
        dcc.Location(id='url', refresh=False),
        
        html.Div(children = [
            
            html.Div(id='page-content')
        
        ])
        
    ])
    
])



page_0_layout = html.Div( children = [
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            html.Img(src='assets/OBP_bas.gif', width=250, height=250),
            
            dmc.Text(
                "iOBPdb Home",
                style={"fontSize": 60},
            ),
            
            html.Img(src='assets/OBP_cart.gif', width=250, height=250)
            
            
        ], style={"marginLeft": 120,"marginBottom":25}
    ),
        
    dmc.Center(children=[
    
        dmc.Grid(children=[
            
            dmc.Col(html.Div(
                
                dmc.Paper(
                    children=[
                
                html.Img(id='schema',),
                
                dmc.Pagination(
                    id='schema_pages',
                    page = 1,
                    total=3,
                    grow=True,
                    align="stretch",
                    siblings=1,
                    
                )
                
                ],shadow="xl",withBorder = True), 
                        

                style={"marginLeft": 250}), span=5),
            
            dmc.Col(html.Div(
                
                dmc.Paper(
                    children=[
                        
                        dmc.Text('',id='explainer_text',style={"fontSize": 18})
                              
                        ],shadow="xl",withBorder = True), 
                
                style={"marginLeft": 50,"marginRight": 170,'height': 320}), span=3),
        
        ],gutter="xs", grow=True ,justify= 'center')
    ])
    
])

         
   
@app.callback(
    Output('schema', 'src'),
    Output('schema', 'style'),
    Output('explainer_text','children'),
    Input('schema_pages', 'page'))

def schema(page): 

    if page == 2:
        
        return ("/assets/s1.png" , {'width':850,
                                    'paddingBottom':40,
                                    'paddingTop':40,
                                    "marginLeft": 0},
                "Odorant binding proteins (OBPs) are a diverse group of small, 10-20 kDa, soluble extracellular target binding proteins found in both terrestrial vertebrates and invertebrates. These studies indicate that there is no shared homology between insect and mammalian OBPs. OBPs in mammals are comprised of a beta barrel type structure, whereas the OBPs in insects are a globular structure comprised of alpha helices. Although OBPs are multifaceted in terms of their potential roles in both insects and mammals alike, they are primarily thought to act as odor transporters, solubilizing volatile organic compounds (VOCs) and pheromones from the surrounding air into the aqueous phase of the odor sensory organ, such as the mucus in the nose or sensillum lymph of an antennae. Once the odors from the environment are solubilized by OBPs, they can then be transported to odor sensing neurons which are coated with olfactory receptor (OR) proteins, which can recognize and bind to specific odors, thus signaling an olfactory response.")
    
    if page == 3:
        
        return ("/assets/s2.png" , {'width':850,"marginLeft": 0, 'paddingBottom':20,'paddingTop':20},
                "The classic insect OBP contains 6 highly conserved cystine residues which form 3 disulphide bonds. However, there are insect OBP variants with fewer than 5 cystines which only from 2 disulphide bonds and are aptly named minus-C OBPs. Conversely there are insect OBPs with 8 or more cystines which are termed as plus-C OBPs. A special variant of plus-C OBPs are atypical OBPs which typically are 20-30 amino acids longer compared to regular insect OBPs and contain 10 or more cystines. Atypical OBPs are also sometimes referred to as two domain OBPs or double domain OBPs which refers to a fusion protein consisting of two OBPs. The variation in cystine count, drastically changes the underlying protein folding, survivability in the extracellular environment and definition of ligand binding pocket.")
    
    if page == 1:
        
        return ("/assets/s4.png" , {'width':850,"marginLeft": 0},
                "Odorant binding proteins, OBPs, are a family of small, globular, extra-cellular proteins which assist in solubilizing volatile organic compounds (VOCs) so they can be internalized and transported by an organism. Since their initial discovery in early eighties, thousands of OBPs have been identified through genome sequencing and characterized by fluorescence ligand binding assays. While a given individual OBP has been studied in the context of their role in specific organism, there has not been studies towards the understanding of the comparative structure-function relations of all known OBPs, primarily due to a lack of a centralized database that incorporates the binding affinity with the structure of all OBPs. Using data obtained from 215 functional studies containing 382 unique OBPs from 91 insect species we created a database, dubbed as iOBPdb, that contains OBP binding affinities for a wide range of VOC targets.")



page_1_layout = html.Div([
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            dmc.Text(
                "Search",
                style={"fontSize": 60},
            )
        ], style={"marginLeft": 120}
    ),
    
    dmc.Paper(children = [
        
        dmc.Center(children = [html.Div([
         
            dmc.Title("Key word search", order=3, align="center"),
            
            dmc.TextInput(id = 'key_search', style={"width": 600,'height':70,"marginTop": 10}),
    
        
        ])], style={"marginLeft": 0,"marginTop": 20}),
    
    
    
        dmc.Center(children = [
            
            dmc.Grid(
                gutter=100, grow=True,
                children=[
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("Compound name search", order=6, align="center"),
                        
                        dmc.Select(id = 'cd_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = unique_list(compact_comp_df,'Compound name'),
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=2),
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("Compound CAS search", order=6, align="center"),
                        
                        dmc.Select(id = 'cas_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = unique_list(ci_df,'CAS-number'),
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=2),
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("OBP name search", order=6, align="center"),
                        
                        dmc.Select(id = 'obp_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = unique_list(oi_df,'Binding Protein Name'),
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=2),
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("OBP Accession search", order=6, align="center"),
                        
                        dmc.Select(id = 'as_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = (unique_list(oi_df,'Accession number / ID') + 
                                           unique_list(oi_df,'UniProtID')), 
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=2)
                ]
            ),
        
        ]),
        
        dmc.Center(children = [
            
            dmc.Grid(
                gutter=40, grow=True,
                children=[
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("Functional group filter", order=6, align="center"),
                        
                        dmc.Select(id = 'func_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = sorted(fg_list),
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=1),
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("Human presence filter", order=6, align="center"),
                        
                        dmc.Select(id = 'hum_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = sorted(hp_list),
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=1),
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("Disease presence filter", order=6, align="center"),
                        
                        dmc.Select(id = 'dis_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = sorted(ds_list),
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=1),
                    
                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("OBP type filter", order=6, align="center"),
                        
                        dmc.Select(id = 'tyoe_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = unique_list(oi_df,'Binding Protein Type'), 
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=1),

                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("OBP species filter", order=6, align="center"),
                        
                        dmc.Select(id = 'spec_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = unique_list(oi_df,'Species'),
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=1),

                    dmc.Col(
                    
                    html.Div([
                     
                        dmc.Title("OBP cystine count filter", order=6, align="center"),
                        
                        dmc.Select(id = 'cys_search', 
                                   searchable=True,
                                   clearable=True,
                                   nothingFound="Not in list",
                                   data = cs_list,
                                   style={"width": 200,'height':70,"marginTop": 10}),
                        
                        ]), span=1),
                    
                ]
            ),
        
        ]),
        
        
    
    ],shadow="xl",withBorder = True, style={"marginLeft": 120,"marginTop": 0}),
    
    
    
    dmc.Center(children = [
        
        dmc.Grid(
            gutter=25, grow=True,
            children=[
                
                dmc.Col(
                
                html.Div([
                    
                    dash_table.DataTable(
                        id={
                            'type': 'data_table_comp',
                            'index': 'comp-table_results'
                        },
                        columns=[{'id': c, 'name': c} for c in compact_comp_df.columns],
                        fill_width=True,
                        editable=True,
                        page_action="native",
                        page_current=0,
                        page_size=5,
                        cell_selectable = True,
                        active_cell=None,
                        style_header={
                        'backgroundColor': 'rgb(30, 30, 30)',
                        'color': 'white'
                        },
                        style_data={
                            'backgroundColor': 'rgb(50, 50, 50)',
                            'color': 'white',"height": 'auto'
                        },
                        style_data_conditional=[                
                            {
                                "if": {"state": "selected"},
                                "backgroundColor": "rgba(50, 50, 50, 1)",
                                "border": "1px solid white"
                            },
                            {'if': {'column_id': 'CAS-number',},
                                'color': '#4dabf7'}
                        ],
                        css=[
                        {"selector": ".dash-spreadsheet-container table", 
                         "rule": '--text-color: #7a4df7 !important'},
                        {'selector': 'td.cell--selected *, td.focused *',
                         'rule': 'text-align: left !important;'
                         }],
                        style_cell={
                            "textAlign": "left",
                            "overflow": "hidden",
                            "textOverflow": "ellipsis",
                            'maxWidth': 0
                        },
                        style_cell_conditional=[
                        {'if': {'column_id': 'id',},
                            'display': 'None',}],
                        style_table={
                            "width": 1200,
                        },
                    
                    )
                    
                    ], style={"paddingLeft":120}), span=5),
                
                dmc.Col(
                
                html.Div([
                    
                    dash_table.DataTable(
                        id={
                            'type': 'data_table_obp',
                            'index': 'obp-table_results'
                        },
                        columns=[{'id': c, 'name': c} for c in compact_OBP_df.columns],
                        fill_width=True,
                        editable=True,
                        page_action="native",
                        page_current=0,
                        page_size=5,
                        cell_selectable = True,
                        active_cell=None,
                        style_header={
                        'backgroundColor': 'rgb(30, 30, 30)',
                        'color': 'white'
                        },
                        style_data={
                            'backgroundColor': 'rgb(50, 50, 50)',
                            'color': 'white',"height": 'auto'
                        },
                        style_data_conditional=[                
                            {
                                "if": {"state": "selected"},
                                "backgroundColor": "rgba(50, 50, 50, 1)",
                                "border": "1px solid white"
                            },
                            {'if': {'column_id': 'Binding Protein Name',},
                                'color': '#4dabf7'}
                        ],
                        css=[
                        {"selector": ".dash-spreadsheet-container table", 
                         "rule": '--text-color: #7a4df7 !important'},
                        {'selector': 'td.cell--selected *, td.focused *',
                         'rule': 'text-align: left !important;'
                         }],
                        style_cell={
                            "textAlign": "left",
                            "overflow": "hidden",
                            "textOverflow": "ellipsis",
                            'maxWidth': 0
                        },
                        style_cell_conditional=[
                        {'if': {'column_id': 'id',},
                            'display': 'None',}],
                        style_table={
                            "width": 1200,
                        },
                    
                    )
                    
                    ], style={"paddingLeft":120}), span=5),
                
            ]
        ),
    
    ], style={"marginLeft": 220,"marginRight": 40}),
        

])

@app.callback(
    Output({'type': 'data_table_comp','index':'comp-table_results'}, 'data'),
    Output({'type': 'data_table_obp','index':'obp-table_results'}, 'data'),
    Input('key_search','value'),
    Input('cd_search','value'),
    Input('cas_search','value'),
    Input('obp_search','value'),
    Input('as_search','value'),
    Input('func_search','value'),
    Input('hum_search','value'),
    Input('dis_search','value'),
    Input('tyoe_search','value'),
    Input('spec_search','value'),
    Input('cys_search','value'))

def table_results(key,cd,cas,obp,asc,func,hum,dis,typ,spec,cys):
 
    s1 = [] ; s2 = []
       
    try:
        
        new_df_comp = compact_comp_df[:-1].loc[compact_comp_df['Compound name'].str.contains(key, case=False)]
        s1 += new_df_comp.to_dict('records')
        
    except:
        
        pass
    
    try:
    
        new_df_comp = compact_comp_df[:-1].loc[compact_comp_df['CAS-number'].str.contains(key, case=False)]
        s1 += new_df_comp.to_dict('records')
        
    except:
        
        pass
    
    try:
    
        cas_list = (ci_df.loc[ci_df[key] == 'Yes'])['CAS-number'].tolist()
        new_df_comp = compact_comp_df[:-1][compact_comp_df['CAS-number'].isin(cas_list)]
        s1 += new_df_comp.to_dict('records')
        
    except:
        
        pass
    
    try:
    
        hum_list = (ci_df.loc[ci_df[key] == 'Yes'])['CAS-number'].tolist()
        new_df_comp = compact_comp_df[:-1][compact_comp_df['CAS-number'].isin(hum_list)]
        s1 += new_df_comp.to_dict('records')
    
    except:
        
        pass
    
    try:
    
        dis_list = (ci_df.loc[ci_df[key] == 'Yes'])['CAS-number'].tolist()
        new_df_comp = compact_comp_df[:-1][compact_comp_df['CAS-number'].isin(dis_list)]
        s1 += new_df_comp.to_dict('records')
    
    except:
        
        pass
    
    try:
    
        new_df_obp =compact_OBP_df[:-1].loc[compact_OBP_df['Binding Protein Name'].str.contains(key, case=False)]
        s2 += new_df_obp.to_dict('records')
    
    except:
        
        pass
    
    try:
    
        new_df_obp =compact_OBP_df[:-1].loc[compact_OBP_df['Accession number / ID'].str.contains(key, case=False)]
        s2 += new_df_obp.to_dict('records')
    
    except:
        
        pass
    
    try:
    
        typ_list = (oi_df.loc[oi_df['Binding Protein Type'] == key])['Binding Protein Name'].tolist()
        new_df_obp = compact_OBP_df[:-1][compact_OBP_df['Binding Protein Name'].isin(typ_list)]
        s2 += new_df_obp.to_dict('records')
    
    except:
        
        pass
    
    try:
    
        spec_list = (oi_df.loc[oi_df['Species'] == key])['Binding Protein Name'].tolist()
        new_df_obp = compact_OBP_df[:-1][compact_OBP_df['Binding Protein Name'].isin(spec_list)]
        s2 += new_df_obp.to_dict('records')
    
    except:
        
        pass
    
    try:
    
        cys_list = (oi_df.loc[oi_df['Cystine count'] == key])['Binding Protein Name'].tolist()
        new_df_obp = compact_OBP_df[:-1][compact_OBP_df['Binding Protein Name'].isin(cys_list)]
        s2 += new_df_obp.to_dict('records')
        
    except:
        
        pass
    
    if cd != None:
        
        new_df_comp = compact_comp_df[:-1].loc[compact_comp_df['Compound name'] == cd]
        s1 += new_df_comp.to_dict('records')
    
    if cas != None:
         
        new_df_comp = compact_comp_df[:-1].loc[compact_comp_df['CAS-number'] == cas]
        s1 += new_df_comp.to_dict('records')
        
    if func != None:
        
        cas_list = (ci_df.loc[ci_df[func] == 'Yes'])['CAS-number'].tolist()
        new_df_comp = compact_comp_df[:-1][compact_comp_df['CAS-number'].isin(cas_list)]
        s1 += new_df_comp.to_dict('records')
        
    if hum != None:
         
        hum_list = (ci_df.loc[ci_df[hum] == 'Yes'])['CAS-number'].tolist()
        new_df_comp = compact_comp_df[:-1][compact_comp_df['CAS-number'].isin(hum_list)]
        s1 += new_df_comp.to_dict('records')
        
    if dis != None:
        
        dis_list = (ci_df.loc[ci_df[dis] == 'Yes'])['CAS-number'].tolist()
        new_df_comp = compact_comp_df[:-1][compact_comp_df['CAS-number'].isin(dis_list)]
        s1 += new_df_comp.to_dict('records')
        
    if obp != None:
        
        new_df_obp =compact_OBP_df[:-1].loc[compact_OBP_df['Binding Protein Name'] == obp]
        s2 += new_df_obp.to_dict('records')
        
    if asc != None:
        
        new_df_obp =compact_OBP_df[:-1].loc[compact_OBP_df['Accession number / ID'] == asc]
        s2 += new_df_obp.to_dict('records')
        
        if s2 == []:
        
            new_df_obp = compact_OBP_df[:-1].loc[compact_OBP_df['UniProtID'] == asc]
            s2 += new_df_obp.to_dict('records')
            
    if typ != None:
        
        typ_list = (oi_df.loc[oi_df['Binding Protein Type'] == typ])['Binding Protein Name'].tolist()
        new_df_obp = compact_OBP_df[:-1][compact_OBP_df['Binding Protein Name'].isin(typ_list)]
        s2 += new_df_obp.to_dict('records')
    
    if spec != None:
        
        spec_list = (oi_df.loc[oi_df['Species'] == spec])['Binding Protein Name'].tolist()
        new_df_obp = compact_OBP_df[:-1][compact_OBP_df['Binding Protein Name'].isin(spec_list)]
        s2 += new_df_obp.to_dict('records')
        
    if cys != None:
        
        cys_list = (oi_df.loc[oi_df['Cystine count'] == cys])['Binding Protein Name'].tolist()
        new_df_obp = compact_OBP_df[:-1][compact_OBP_df['Binding Protein Name'].isin(cys_list)]
        s2 += new_df_obp.to_dict('records')
    
    return s1, s2





page_2_layout = html.Div(children =[
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            dmc.Text(
                "Browse Compounds",
                style={"fontSize": 60},
            )
        ], style={"marginLeft": 120}
    ),
    
    dmc.Center(children=[
    
        dash_table.DataTable(
            id={
                'type': 'data_table_comp',
                'index': 'comp-table'
            },
            columns=[{'id': c, 'name': c} for c in compact_comp_df.columns],
            data = compact_comp_df.to_dict('records'),
            editable=True,
            page_action="native",
            page_current=0,
            active_cell=None,
            page_size=25,
            cell_selectable = True,
            style_header={
            'backgroundColor': 'rgb(30, 30, 30)',
            'color': 'white'
            },
            style_data={
                'backgroundColor': 'rgb(50, 50, 50)',
                'color': 'white',"height": 'auto'
            },
            style_data_conditional=[                
                {
                    "if": {"state": "selected"},
                    "backgroundColor": "rgba(50, 50, 50, 1)",
                    "border": "1px solid white"
                },
                {'if': {'column_id': 'CAS-number',},
                    'color': '#4dabf7'}
            ],
            css=[
            {"selector": ".dash-spreadsheet-container table", 
             "rule": '--text-color: #7a4df7 !important'},
            {'selector': 'td.cell--selected *, td.focused *',
             'rule': 'text-align: left !important;'
             }],
            style_cell={
                "textAlign": "left",
                "overflow": "hidden",
                "textOverflow": "ellipsis",
                'maxWidth': 0
            },
            style_cell_conditional=[
            {'if': {'column_id': 'id',},
                'display': 'None',}],
            style_table={
                "width": 1200,
            },
        
        ),
        
    ], style={"paddingLeft":120})

])



page_3_layout = html.Div(children =[
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            dmc.Text(
                "Browse OBPs",
                style={"fontSize": 60},
            )
        ], style={"marginLeft": 120}
    ),
    
    dmc.Center(children=[
    
        dash_table.DataTable(
            id={
                'type': 'data_table_obp',
                'index': 'obp-table'
            },
            columns=[{'id': c, 'name': c} for c in compact_OBP_df.columns],
            data = compact_OBP_df.to_dict('records'),
            fill_width=True,
            editable=True,
            page_action="native",
            page_current=0,
            page_size=25,
            cell_selectable = True,
            active_cell=None,
            style_header={
            'backgroundColor': 'rgb(30, 30, 30)',
            'color': 'white'
            },
            style_data={
                'backgroundColor': 'rgb(50, 50, 50)',
                'color': 'white',"height": 'auto'
            },
            style_data_conditional=[                
                {
                    "if": {"state": "selected"},
                    "backgroundColor": "rgba(50, 50, 50, 1)",
                    "border": "1px solid white"
                },
                {'if': {'column_id': 'Binding Protein Name',},
                    'color': '#4dabf7'}
            ],
            css=[
            {"selector": ".dash-spreadsheet-container table", 
             "rule": '--text-color: #7a4df7 !important'},
            {'selector': 'td.cell--selected *, td.focused *',
             'rule': 'text-align: left !important;'
             }],
            style_cell={
                "textAlign": "left",
                'width': '{}%'.format(len(compact_OBP_df.columns)),
                'textOverflow': 'ellipsis',
                'overflow': 'hidden'
            },
            style_cell_conditional=[
            {'if': {'column_id': 'id',},
                'display': 'None',}],
            style_table={
                "width": 1200,
            },
        
        ),
        
    ], style={"paddingLeft":920,"width": 900})

])



df = pd.DataFrame(columns=['Compound Name',
                           'VOC (CAS#)*',
                           'Compound SMILES',
                           'Binding Affinity (Kd)*'], index=range(5))

page_4_layout = html.Div([
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            dmc.Text(
                "Data Entry Form",
                style={"fontSize": 60},
            )
        ], style={"marginLeft": 120,"marginBottom": 32}
    ),
    
    dmc.Grid(
        children=[
            dmc.Col(html.Div(dmc.Stack(children=[
                
                dmc.Text(
                    "OBP Data:",
                    style={"fontSize": 30},
                ),
                
                dmc.TextInput(label="Submitter:", 
                              style={"width": 400}, 
                              placeholder="Your Name"),
                
                dmc.TextInput(label="Email*:", 
                              style={"width": 400}, 
                              placeholder="Your Email", 
                              ), #error="Enter a valid email"
                
                dmc.TextInput(label="DOI*:", 
                              style={"width": 400}, 
                              placeholder="Publication DOI"),
                
                dmc.TextInput(label="Protein name*:", 
                              style={"width": 400}, 
                              placeholder="Name of OBP, CSP, PBP etc."),
                
                dmc.TextInput(label="Species name:", 
                              style={"width": 400}, 
                              placeholder="Name of species"),
                
                dmc.TextInput(label="Online deposition of protein (UniProt, GenBank, RCSB etc.):", 
                              style={"width": 400}, 
                              placeholder="IDs, URL, accesion numbers"),
                
                dmc.Textarea(
                    label="Protein sequence*:",
                    placeholder="Amino acid sequence for full protein",
                    style={"width": 400},
                    autosize=True,
                    minRows=3,
                ),
                
                dmc.Textarea(
                    label="Additional comments:",
                    placeholder="Description of protein, type of binding experiment, other notes",
                    style={"width": 400},
                    autosize=True,
                    minRows=3,
                ),
                
            ], style={"marginLeft": 240})), span=5),
            
            dmc.Col(html.Div([
                
                dmc.Text(
                    "OBP-VOC Binding Data:",
                    style={"fontSize": 30,"marginBottom": 32},
                ),
                
                dash_table.DataTable(
                    id='data-table',
                    columns=[{'name': col, 'id': col} for col in df.columns],
                    data=df.to_dict('records'),
                    editable=True,
                    row_deletable=True,
                    page_size=15,
                    style_header={
                    'backgroundColor': 'rgb(30, 30, 30)',
                    'color': 'white'
                    },
                    style_data={
                    'backgroundColor': 'rgb(50, 50, 50)',
                    'color': 'white',
                    'width': '155px'
                    },
                    style_cell={'textAlign': 'left',
                                'overflow': 'hidden',        
                                'textOverflow': 'ellipsis',
                                'maxWidth': '155px',
                                },
                    
                    style_data_conditional=[
                        {
                            # Style for editable cells
                            'if': {'column_editable': True},
                            'color': 'white',
                            'backgroundColor': 'rgb(50, 50, 50)',
                        },
                        {
                            # Style for selected cells
                            'if': {'state': 'selected'},
                            'color': 'white',
                            'backgroundColor': 'rgba(255, 255, 255, 1)',
                        },
                        {
                            # Style for active (currently editing) cells
                            'if': {'state': 'active'},
                            'color': 'white',
                            'backgroundColor': 'rgba(255, 255, 255, 1)',
                        },
                        {'if': {'column_id': 'Compound name'},
                             'width': '25%'},
                        {'if': {'column_id': 'VOC (CAS#)*'},
                             'width': '25%'},
                        {'if': {'column_id': 'Compound SMILES'},
                             'width': '25%'},
                        {'if': {'column_id': 'Binding Affinity (kD)*'},
                             'width': '25%'},
                    ]),
                
                html.Button('Add Row', id='add-row-button', n_clicks=0),
                html.Button('Submit', id='submit-button', n_clicks=0),
                html.Div(id='output')
                
            ], style={"marginRight": 60}), span = 7),
        ],gutter="xl")

])

@app.callback(
    dash.dependencies.Output('data-table', 'data'),
    [dash.dependencies.Input('add-row-button', 'n_clicks')],
    [dash.dependencies.State('data-table', 'data'),
     dash.dependencies.State('data-table', 'columns')]
)
def add_row(n_clicks, table_data, table_columns):
    if n_clicks > 0:
        table_data.append({c['id']: '' for c in table_columns})
    return table_data

@app.callback(
    dash.dependencies.Output('output', 'children'),
    [dash.dependencies.Input('submit-button', 'n_clicks')]
)
def update_output(n_clicks):
    if n_clicks > 0:
        return 'Thank you for submitting your information'




files = {
    'OBP_info.csv': {'path': 'OBP_info.csv', 
               'text': 'The characteristics and information pertinent to every odorant binding protein (OBP) are stored in a data table on OBP_info.csv'},
    'Compound_info.csv': {'path': 'Compound_info.csv', 
               'text': 'The specific attributes of every volatile organic compound (VOC) are found in a data table stored on Compound_info.csv'},
    'Compound_OBP_binding.csv': {'path': 'Compound_OBP_binding.csv', 
               'text': 'The binding affinity of every recorded OBP-VOC interaction is found in a data table stored on Compound_OBP_binding.csv'},
    'iOBPdb_PDBs_alphafold_denovo.zip': {'path': 'iOBPdb_PDBs_alphafold_denovo.zip', 
               'text': 'All de novo generated PDB structures with AlphaFold V2.0 is stored in iOBPdb_PDBs_alphafold_denovo.zip'},
    'iOBPdb_full_site.zip': {'path': 'iOBPdb_full_site.zip', 
               'text': 'All previous files stored in a unified zip file'}}

page_5_layout = html.Div([
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            dmc.Text(
                "Download Data Files",
                style={"fontSize": 60},
            ),
            
            
        ], style={"marginLeft": 120}
    ),
    
    html.Div(style={'height': '30px'}),
    
    dmc.Center(
        
        children=[
    
            html.Div([
                
                html.Table([
                    html.Tr([
                        html.Td([
                            html.Div(
                                files[file_name]['text'],
                                style={'font-size': '24px'}
                            ),
                            
                            html.Div([
                            
                                html.H1('Download:', style={'font-size': '24px',
                                                           'display': 'inline-block'}),
                                
                                html.A(
                                file_name,
                                href=files[file_name]['path'],
                                download=file_name,
                                style={'color': '#4dabf7',
                                       'font-size': '24px'})
                                
                                ]),
                            
                            html.Div(style={'height': '30px'})
                        ],style={'margin': 'auto'})
                    ]) for file_name in files ],
                    
                style={'margin': 'auto'})
        
        ], style={"marginLeft": 480,
                  "marginRight": 340})])
])






page_6_layout = html.Div([
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            dmc.Text(
                "Compound Info",
                id='Compound_Info',
                style={"fontSize": 60},
            )
            
        ], style={"marginLeft": 120}
        
    ),

    dmc.Grid(children=[
        
        dmc.Col(html.Div([
            
            html.Img(id='molimg', 
                     style={'height':'300 px', 
                            'width':'300 px',
                            'margin':0})
            
            ]), span=1),
        
        dmc.Col(html.Div(
            
            html.Div([
            
                dmc.Title("Compound alternate name(s)", order=6),
                
                dmc.Text(
                    "-",
                    id='molname',
                    style={"fontSize": 16},
                ),
                
                dmc.Title("Compound CAS", order=6),
                
                dmc.Text(
                    "-",
                    id='molcas',
                    style={"fontSize": 16},
                ),
                
                dmc.Title("Compound formula", order=6),
                
                dmc.Text(
                    "-",
                    id='molform',
                    style={"fontSize": 16},
                ),
                
                dmc.Title("Compound molecular weight", order=6),
                
                dmc.Text(
                    "-",
                    id='molmass',
                    style={"fontSize": 16},
                ),
                
                dmc.Title("Compound smiles", order=6),
                
                dmc.Text(
                    "-",
                    id='molsmile',
                    style={"fontSize": 16},
                ),
                
            ], style={"marginBottom": 15})), span=3),
        
        dmc.Col(html.Div(
            
            html.Div([
            
                dmc.Title("Compound functional groups", order=6),
                
                dmc.Text(
                    "-",
                    id='molfunc',
                    style={"fontSize": 16},
                ),
                
                dmc.Title("Compound found in human", order=6),
                
                dmc.Text(
                    "-",
                    id='molhum',
                    style={"fontSize": 16},
                ),
                
                dmc.Title("Compound found in disease", order=6),
                
                dmc.Text(
                    "-",
                    id='moldis',
                    style={"fontSize": 16},
                ),
                
            ])), span=3),
            

    ],gutter="xs", grow=True ,justify= 'center', style={"marginLeft": 320,"marginRight": 230}),
    
    
    
    dmc.Center(children=[
    
        dash_table.DataTable(
            id={
                'type': 'data_table_obp',
                'index': 'comp_spec-table'
            },
            columns=[{'name':'initial','id':'initial'}],
            fill_width=True,
            editable=True,
            page_action="native",
            page_current=0,
            page_size=15,
            cell_selectable = True,
            active_cell=None,
            style_header={
            'backgroundColor': 'rgb(30, 30, 30)',
            'color': 'white'
            },
            style_data={
                'backgroundColor': 'rgb(50, 50, 50)',
                'color': 'white',"height": 'auto'
            },
            style_data_conditional=[                
                {
                    "if": {"state": "selected"},
                    "backgroundColor": "rgba(50, 50, 50, 1)",
                    "border": "1px solid white"
                },
                {'if': {'column_id': 'Binding Protein Name',},
                    'color': '#4dabf7'}
            ],
            css=[
            {"selector": ".dash-spreadsheet-container table", 
             "rule": '--text-color: #7a4df7 !important'},
            {'selector': 'td.cell--selected *, td.focused *',
             'rule': 'text-align: left !important;'
             }],
            style_cell={
                "textAlign": "left",
                "overflow": "hidden",
                "textOverflow": "ellipsis",
                'maxWidth': 0
            },
            style_cell_conditional=[
            {'if': {'column_id': 'id',},
                'display': 'None',}],
            style_table={
                "width": 1200,
            },
        
        ),
        
    ], style={"paddingLeft":920,"width": 900})

])

@app.callback(
    Output('Compound_Info', 'children'),
    Output('molimg', 'src'),
    Output('molname', 'children'),
    Output('molcas', 'children'),
    Output('molform', 'children'),
    Output('molmass', 'children'),
    Output('molsmile', 'children'),
    Output('molfunc','children'),
    Output('molhum','children'),
    Output('moldis','children'),
    Input('row_selected_comp', 'data')
)

def update_img(cas):

    row = ci_df[ci_df['CAS-number'] == cas].index[0]
    smile = ci_df[:-1].at[row, "Smiles"]
    name = ci_df.at[row, "Compound name"]
     
    title_name = name.split(' /', 1)[0]
    try:
        title_name = title_name.capitalize()
    except:
        pass
    alt_name = name.split(' /', 1)[1:]
    if alt_name == []:
        alt_name = '-'
    else:
        try:
            alt_name = alt_name.capitalize()
        except:
            pass
    title = title_name + ' Info'

    form = ci_df.at[row, "Molecular formula"]
    mass = str(float(Formula(form).mass)) + ' g/mol'
    temp_url = 'https://cactus.nci.nih.gov/chemical/structure/' + str(cas) + '/image?format=png&bgcolor=transparent&antialiasing=0&atomcolor=white&bondcolor=white&width=250&height=250&linewidth=2&symbolfontsize=16'
    
    try:
        
        urllib.request.urlretrieve(temp_url)
    
    except:
    
        temp_url = 'https://cactus.nci.nih.gov/chemical/structure/' + str(name) + '/image?format=png&bgcolor=transparent&antialiasing=0&atomcolor=white&bondcolor=white&width=250&height=250&linewidth=2&symbolfontsize=16'
        
        try:
        
            urllib.request.urlretrieve(temp_url)
            
        except:
    
            temp_url = 'https://cactus.nci.nih.gov/chemical/structure/' + str(smile) + '/image?format=png&bgcolor=transparent&antialiasing=0&atomcolor=white&bondcolor=white&width=250&height=250&linewidth=2&symbolfontsize=16'  
        
    new_df = ci_df[:-1].loc[ci_df['CAS-number'] == cas]
    
    func_list = new_df.iloc[:, 5:-24]
    func_list = func_list.dropna(axis=1, how='all')
    func_list = str(list(func_list))
    if func_list != []:
        func_list = func_list[2:-2].replace("', '", ' / ')
    else:
        func_list ='-'
    
    hum_list = new_df.iloc[:, -24:-17]
    hum_list = hum_list.dropna(axis=1, how='all')
    hum_list = str(list(hum_list))
    if hum_list != '[]':
        hum_list = hum_list[2:-2].replace("', '", ' / ')
    else:
        hum_list ='-'
    
    dis_list = new_df.iloc[:, -17:]
    dis_list = dis_list.dropna(axis=1, how='all')
    dis_list = str(list(dis_list))
    if dis_list != '[]':
        dis_list = dis_list[2:-2].replace("', '", ' / ')
    else:
        dis_list ='-'
    
    try:
        
        return title, temp_url, alt_name, cas, form, mass, smile, func_list, hum_list,dis_list
        
    except:
        
        pass
        
     

@app.callback(Output({'type': 'data_table_obp','index':'comp_spec-table'}, 'data'),
              Output({'type': 'data_table_obp','index':'comp_spec-table'}, 'columns'),
              Input('row_selected_comp', 'data'))

def display_table_comp(cas):
    
    new_df = cob_df[:-1].loc[cob_df['CAS-number'] == cas] 
    new_df = new_df.dropna(axis=1, how='all').T
    new_df.reset_index(inplace=True)
    new_df = new_df[2:]  
    new_df.columns = ['Binding Protein Name', 'Binding affinity (Ki uM)']
    sources_df = oi_df[['Binding Protein Name', 'Source(s)']]
    new_df = pd.merge(new_df, sources_df)
    new_df['id'] = new_df['Binding Protein Name']
    columns=[{'id': c, 'name': c} for c in new_df.columns]
    data = new_df.to_dict('records')
    
    return data, columns

page_7_layout = html.Div([
    
    Navbar,
    
    dmc.Center(
        
        children=[
            
            dmc.Text(
                "OBP Info",
                id = 'OBP_info',
                style={"fontSize": 60},
            )
        ], style={"marginLeft": 120}
    ),
    
    dmc.Center(children=[
    
        dmc.Grid(children=[
            
            dmc.Col(html.Div(
                
                html.Div([
                
                    dmc.Title("OBP PDBs", order=6,style={"paddingTop": 90}),
                    
                    dmc.Select(
                        id="dropdown-pdb-OBP",
                        value = os.path.join(DATAPATH, '1N8U.pdb'),
                        data =pdb_options, style={"width": 130, "marginBottom": 20},
                        
                    ),
                    
                    dmc.Title("PDB style", order=6),
                    
                    dmc.Select(
                        id="dropdown-pdb-style",
                        value="Cartoon",
                        data=[
                            {"value": "Cartoon", "label": "Cartoon"},
                            {"value": "BAS", "label": "Ball and Stick"},
                        ], style={"width": 130, "marginBottom": 20},
                        
                    ),
                    
                    dmc.Title("Color Select", order=6),
                    
                    dbc.Input(
                        type="color",
                        id="pdb_color",
                        value="#ffffff",
                        style={"width": 130,
                               "height": 35,
                               "marginBottom": 20}),
                    
                    dmc.Title("PDB Source", order=6),
                    
                    html.Div(
                        dcc.Link('-', id='pdbsource', href='',
                                 target='_blank',
                                 style={'textDecoration': 'underline',
                                        'color': '#4dabf7',
                                        'cursor': 'pointer'}),
                        style={'width': '130px'}
                    )
                    
                ]), style={"marginLeft": 80,"marginRight": 25}), span=1),
                
                
            dmc.Col(html.Div(
                
                dmc.Skeleton(visible=True,
                             animate=True,
                             height=500,
                             width = 500,
                             id='load_pdb',
                             children=html.Div(
                                id='mol3d-biomolecule-OBP',
                                children=[]
                                )), style={"marginLeft": 0,
                                           "marginRight": 25,
                                           "marginTop": 10}
            ), span=3),
            
            dmc.Col(html.Div(
                
                html.Div([
        
                    dmc.Title("OBP full sequence", order=6),
                    
                    dashbio.SequenceViewer(
                        id='OBPseq',
                        showLineNumbers=False,
                        charsPerLine=50,),
                    
                    dmc.Title("OBP cleaved sequence ", order=6),
                    
                    dashbio.SequenceViewer(
                        id='OBPseq_sp',
                        showLineNumbers=False,
                        charsPerLine=50,),
        
                    dmc.Grid(children=[ 
                        
                        dmc.Col(html.Div([
                        
                            dmc.Title("Source paper", order=6),
                        
                            dmc.Text(
                                "-",
                                id='obpsource',
                                style={"fontSize": 16}),
                        
                            dmc.Title("GenBank Accession", order=6),
                            
                            dmc.Text(
                                "-",
                                id='obpac',
                                style={"fontSize": 16},
                            ),
                            
                            dmc.Title("UniProt Accession", order=6),
                            
                            dmc.Text(
                                "-",
                                id='obpup',
                                style={"fontSize": 16},
                            ),
                            
                            dmc.Title("OBP species", order=6),
                        
                            dmc.Text(
                                "-",
                                id='obpspec',
                                style={"fontSize": 16}),
                            
                            ], style={"marginRight": 0}), span=4),
                        
                        dmc.Col(html.Div([
                            
                            dmc.Title("OBP type", order=5),
                            
                            dmc.Text(
                                "-",
                                id='obptype',
                                style={"fontSize": 16},
                            ),
                            
                            dmc.Title("Molecular Weight", order=6),
                            
                            dmc.Text(
                                "-",
                                id='obpcys',
                                style={"fontSize": 16},
                            ),
                            
                            dmc.Title("OBP PDB entries", order=6),
                            
                            dmc.Text(
                                "-",
                                id='obppdbrcsb',
                                style={"fontSize": 16}),
                            
                            dmc.Text(
                                "-",
                                id='obppdbsm',
                                style={"fontSize": 16}),
                            
                            dmc.Text(
                                "-",
                                id='obppdbalpha',
                                style={"fontSize": 16}),
                            
                            
                            ], style={"marginLeft": 80}), span=7),
                        
                        
                        ])
                ],style={'width':520})
            
            ), span=4)
            
        ],gutter="xs", grow=True ,justify= 'center')
            
    ]),
    
    dmc.Center(children=[
    
        dash_table.DataTable(
            id={
                'type': 'data_table_comp',
                'index': 'obp_spec-table'
            },
            fill_width=True,
            editable=True,
            page_action="native",
            page_current=0,
            page_size=8,
            cell_selectable = True,
            active_cell=None,
            style_header={
            'backgroundColor': 'rgb(30, 30, 30)',
            'color': 'white'
            },
            style_data={
                'backgroundColor': 'rgb(50, 50, 50)',
                'color': 'white',"height": 'auto'
            },
            style_data_conditional=[                
                {
                    "if": {"state": "selected"},
                    "backgroundColor": "rgba(50, 50, 50, 1)",
                    "border": "1px solid white"
                },
                {'if': {'column_id': 'CAS-number',},
                    'color': '#4dabf7'}
            ],
            css=[
            {"selector": ".dash-spreadsheet-container table", 
             "rule": '--text-color: #7a4df7 !important'},
            {'selector': 'td.cell--selected *, td.focused *',
             'rule': 'text-align: left !important;'
             }],
            style_cell={
                "textAlign": "left",
                "overflow": "hidden",
                "textOverflow": "ellipsis",
                'maxWidth': 0
            },
            style_cell_conditional=[
            {'if': {'column_id': 'id',},
                'display': 'None',}],
            style_table={
                "width": 1200,
                "marginTop": 25
            },
        
        ),
        
    ], style={"paddingLeft":920,"width": 900})

])

@app.callback(
    Output('OBP_info', 'children'),
    Output('OBPseq', 'sequence'),
    Output('OBPseq_sp', 'sequence'),
    Output('OBPseq', 'coverage'),
    Output('OBPseq_sp', 'coverage'),
    Output('obpsource', 'children'),
    Output('obpac', 'children'),
    Output('obpup', 'children'),
    Output('obpspec', 'children'),
    Output('obptype', 'children'),
    Output('obppdbrcsb', 'children'),
    Output('obppdbsm', 'children'),
    Output('obppdbalpha', 'children'),
    Output('obpcys', 'children'),
    Input('row_selected_OBP', 'data')
)

def update_seq(obp):
    
    row = oi_df[oi_df['Binding Protein Name'] == obp].index[0]
    obp_title = oi_df[:-1].at[row, "Binding Protein Name"] + ' Info'
    obp_source = oi_df[:-1].at[row, "Source(s)"]
    
    if len(obp_source) > 32:
        
        obp_source = obp_source[:30] + ' ' + obp_source[30:]
        
        if len(obp_source) > 62:
            
            obp_source = obp_source[:60] + ' ' + obp_source[60:]
            
            if len(obp_source) > 92:
                
                obp_source = obp_source[:90] + ' ' + obp_source[90:]
    
    obp_ac = oi_df[:-1].at[row, "Accession number / ID"]
    obp_up = oi_df[:-1].at[row, "UniProtID"]
    obp_spec = oi_df[:-1].at[row, "Species"]
    obp_type = oi_df[:-1].at[row, "Binding Protein Type"]
    
    obp_pdbrcsb = 'RCSB: ' + oi_df[:-1].at[row, "RCSB entry"]
    obp_pdbsm = '\nSWISS-MODEL: ' + oi_df[:-1].at[row, "SwissModel PDB"]
    obp_pdbalpha = '\nAlphaFold: ' + oi_df[:-1].at[row, "AlphaFold PDB"]
    
    obp_seq = oi_df[:-1].at[row, "AA Sequence"]
    
    cysteine_indices = [i for i, res in enumerate(obp_seq) if res == 'C']
    obp_seq_cov = [{'start': idx, 'end': idx+1, 
                    'color': 'yellow'} for idx in cysteine_indices]
    non_coverage = [{'start': idx, 'end': idx+1, 
                     'color': 'white'} for idx in range(len(obp_seq)) if idx not in cysteine_indices]

    obp_seq_cov += non_coverage
    
    obp_seq_sp = oi_df[:-1].at[row, "AA Sequence W/O signal peptide"]
    cysteine_indices = [i for i, res in enumerate(obp_seq_sp) if res == 'C']
    obp_seq_sp_cov = [{'start': idx, 'end': idx+1, 
                       'color': 'yellow'} for idx in cysteine_indices]
    non_coverage = [{'start': idx, 'end': idx+1, 
                     'color': 'white'} for idx in range(len(obp_seq_sp)) if idx not in cysteine_indices]

    obp_seq_sp_cov += non_coverage
    
    protein_analysis = ProteinAnalysis(obp_seq_sp)
    obp_cys = protein_analysis.molecular_weight() / 1000
    obp_cys = '{:.3f}'.format(obp_cys) + ' kDa'
    
    try:
         
        return (obp_title, obp_seq, obp_seq_sp, obp_seq_cov, obp_seq_sp_cov,
                obp_source, obp_ac,obp_up,obp_spec,
                obp_type,obp_pdbrcsb,obp_pdbsm,obp_pdbalpha,obp_cys)
        
    except:
        
        pass
    
@app.callback(Output({'type': 'data_table_comp','index':'obp_spec-table'}, 'data'),
              Output({'type': 'data_table_comp','index':'obp_spec-table'}, 'columns'),
              Input('row_selected_OBP', 'data'))

def display_table_cobp(obp_name):
    
    new_df = cob_df[['CAS-number','Compound name',obp_name]]
    new_df.columns = [*new_df.columns[:-1], "Binding affinity (Ki uM)"]
    new_df = (new_df.dropna())
    new_df['id'] = new_df['CAS-number'] 
    
    columns=[{'id': c, 'name': c} for c in new_df.columns]
    
    data = new_df.to_dict('records')
    
    return data, columns



@app.callback(Output('page-content', 'children'),
              Output('url', 'pathname'),
              Output('row_selected_comp', 'data'),
              Output('row_selected_OBP', 'data'),
              Input('url', 'pathname'),
              Input('row_selected_comp', 'data'),
              Input('row_selected_OBP', 'data'),
              Input({'type': 'data_table_comp','index':ALL}, 'active_cell'),
              Input({'type': 'data_table_obp','index':ALL}, 'active_cell'))
   
def display_page(pathname,cas,obp,cell0,cell1):
    
    ctx = dash.callback_context 
    ctx_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    try:
        
        if cell0[0]["column_id"] == 'CAS-number':
    
            cas = cell0[0]["row_id"]
        
    except: 
        
        pass 
    
    try:
        
        if cell1[0]["column_id"] == 'Binding Protein Name':
    
            obp = cell1[0]["row_id"]
        
    except: 
        
        pass 
    
    try:
        
        if pathname == '/About':
            
            if ctx_id != 'url':
                
                return dash.no_update 
            
            return page_0_layout,pathname,cas,obp
        
        elif pathname == '/Search':
            
            if ctx_id != 'url':
                
                if cas != None:
                    
                    path = '/' + cas
                     
                    return page_6_layout,path,cas,obp
                
                elif obp != None: 
                    
                    path = '/' + obp
                          
                    return page_7_layout,path,cas,obp
                
                return dash.no_update 
                
            return page_1_layout,pathname,cas,obp
        
        elif pathname == '/BrowseCompounds':
        
            if cas != None:
                
                path = '/' + cas
                 
                return page_6_layout,path,cas,obp
            
            if cell0 == [None]:
                
                return dash.no_update
        
            return page_2_layout,pathname,cas,obp
        
        elif pathname == '/BrowseOBPs':
            
            if obp != None: 
                
                path = '/' + obp
                      
                return page_7_layout,path,cas,obp
            
            if cell1 == [None]:
                
                return dash.no_update
            
            return page_3_layout,pathname,cas,obp
        
        elif pathname == '/'  + str(obp):
            
            if cell0 == [None]:
                
                return dash.no_update
            
            if cas != None:
                
                path = '/' + cas
                
                return page_6_layout,path,cas,obp
        
        elif pathname == '/' + str(cas):
            
            if cell1 == [None]:
                
                return dash.no_update
            
            if obp != None: 
                
                path = '/' + obp
                
                return page_7_layout,path,cas,obp
        
        elif pathname == '/DataEntry':
            
            return page_4_layout,pathname,cas,obp
        
        elif pathname == '/Download':
            
            return page_5_layout,pathname,cas,obp
            
        else:
            
            return page_0_layout,pathname,cas,obp
        
    except:
        
        return dash.no_update 
    
    
    
@app.callback(
    Output('dropdown-pdb-OBP', 'value'),
    Output('dropdown-pdb-OBP', 'data'),
    Input('row_selected_OBP', 'data'),
)

def OBP_dropdown(obp):

    row = oi_df[oi_df['Binding Protein Name'] == obp].index[0]
    DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PDB_iodbp')
    pdb_options = []
    
    obp_pdbrcsb = (oi_df[:-1].at[row, "RCSB entry"]).split(" / ")

    for PDB in obp_pdbrcsb:
        
        if PDB != '-':
        
            pdb_options += [{
                'label': PDB,
                'value': os.path.join(DATAPATH, PDB +'.pdb')
            }]
        
    DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SWM_pdbs')
    obp_pdbsm = oi_df[:-1].at[row, "SwissModel PDB"]
    
    if obp_pdbsm != '-': 
    
        pdb_options += [{
            'label': 'SwissModel',
            'value': os.path.join(DATAPATH, obp_pdbsm)
        }]
    
    DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'alpha_pdbs')
    obp_pdbalpha = oi_df[:-1].at[row, "AlphaFold PDB"]
      
    pdb_options += [{ 
        'label': 'AlphaFold',
        'value': os.path.join(DATAPATH, obp_pdbalpha)
    }]
    value = os.path.join(DATAPATH, obp_pdbalpha)
    
    return value, pdb_options

    
@app.callback(
    Output('mol3d-biomolecule-OBP', 'children'),
    Output('pdbsource', 'children'),
    Output('pdbsource', 'href'),
    Output('load_pdb','visible'),
    Input('dropdown-pdb-OBP', 'value'),
    Input('dropdown-pdb-OBP', 'data'),
    Input('dropdown-pdb-style', 'value'),
    Input("pdb_color", "value"),
    Input('obpup', 'children')
)
    
def use_upload_obp(pdb,options,style,color,uniprotID):

    selected_option = next((option for option in options if option['value'] == pdb), None)

    if pdb is not None:
        
        copy2(pdb, './str.pdb')
        fname = './str.pdb'

    else:
         copy2((os.path.join(DATAPATH, '1N8U.pdb')), './str.pdb')
         fname = './str.pdb'
        
    pdbsource = selected_option["label"]
    
    if pdbsource == None:
        
        text = ''
        link = ''
    
    elif pdbsource == 'AlphaFold':
    
        text = 'PDB generated with AlphaColab'
        link = 'https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb'
        
    elif pdbsource == 'SwissModel':
    
        text = 'PDB for ' + uniprotID + ' found on Swiss Model server'
        #'generated with Swiss Model server'
        link = 'https://swissmodel.expasy.org/repository/uniprot/' + uniprotID
        
    else:
        
        text = 'PDB for ' + pdbsource + ' found on RCSB'
        link = 'https://www.rcsb.org/structure/' + pdbsource
        
    # elif pdbsource == 'SwissModel'
        
    # Create the model data from the decoded contents
    pdb = parser.PdbParser(fname)
    mdata = pdb.mol3d_data()
    
    CHAIN_COLORS = {
    "A": color,
    "B": color,
    "C": color,
    "D": color,
    "E": color,
    "F": color,
    "G": color,
    "H": color,
    "I": color,
    "J": color,
    "K": color,
    "L": color,
    "M": color,
    "N": color,
    "O": color,
    "P": color,
    "R": color,
    "S": color,
}

    if style == 'BAS':

        style = sparser.create_mol3d_style(
                mdata.get("atoms"),
                'stick',
                'atom')
    
    else:
        
        style = sparser.create_mol3d_style(
                mdata.get("atoms"),
                'cartoon',
                'chain',
                CHAIN_COLORS)

    return dash_bio.Molecule3dViewer(
            id='mol-3d',
            selectionType='atom',
            modelData=mdata,
            styles=style,
            selectedAtomIds=[],
            backgroundOpacity='0',
            atomLabelsShown=False,
            zoom = {'factor': 1.35, 'animationDuration': 0, 'fixedPath': False},
            style = { 'height': 500, 'width': 500}), text, link, False





if __name__ == '__main__':
    app.run_server(debug=True)