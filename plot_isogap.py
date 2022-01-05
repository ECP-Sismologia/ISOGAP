# /**
#  * @author Emmanuel Castillo
#  * @email ecastillot@unal.edu.co
#  * @create date 2021-08-19 09:33:05
#  * @modify date 2021-08-19 09:33:05
#  * References:
        # https://stackoverflow.com/questions/36556890/restructuring-pandas-dataframe-into-meshgrid-for-basemap (basemap)
        # https://www.pygmt.org/dev/gallery/images/grdclip.html#sphx-glr-gallery-images-grdclip-py (pygmt)
        # https://forum.generic-mapping-tools.org/t/map-of-seismic-stations-using-pygmt/828 (pygmt to maps)
#  */
import os
import pygmt
import xarray
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import geopandas as gpd
import ast
from itertools import takewhile
pygmt.config(FORMAT_GEO_MAP="ddd.xx")

class ShapeObject():
    def __init__(self,data,label,
                **kwargs):
        " kwargs of pygmt plot"
        self.data = data
        self.label = label
        self.__dict__.update(**kwargs)

def get_shape(shp_path,name_dict=None,
                proj="EPSG:4326",**kwargs):
    gdf = gpd.read_file(shp_path)
    if proj != "EPSG:4326":
        gdf = gdf.to_crs(proj)
    else: 
        gdf = gdf.to_crs("EPSG:4326")

    if name_dict != None: 
        data = []   
        for key,values in name_dict.items():
            for value in values:
                df = gdf[gdf[key] == value]

                print(df)
                data.append(df)
    else: 
        data = [gdf]

    

    return data

def get_csv_comments(csv_path):
    with open(csv_path, 'r') as fobj:

        headiter = takewhile(lambda s: s.startswith(f'#'), fobj)
        headers = list(headiter)

    comments = {}
    for comment in headers:
        comment = comment.strip('#').strip('\n').strip()
        key,value = comment.split("=")
        comments[key] = value

    return comments

def get_stations(csv_path,off_stations=[]):
    df = pd.read_csv(csv_path)

    df["status"] = df["STA"].apply(lambda x:  False if x in off_stations else True)
    on_df = df[ df["status"]==True]
    off_df = df[ df["status"]==False]

    stations = {}
    for status,df in [("ON",on_df),("OFF",off_df)]:
        sta = df["STA"].to_numpy()
        lon = df["LON"].to_numpy()
        lat = df["LAT"].to_numpy()
        stations[status] = [sta,lon,lat]
    return stations

def get_meshgrid(csv_path):
    comments= get_csv_comments(csv_path)
    print(comments)
    if not comments:
        off_stations = []
    else:
        off_stations = ast.literal_eval(comments["off_stations"])

    df = pd.read_csv(csv_path, comment='#')
    pivot_df = df.pivot(index='LAT', columns='LON', values='GAP')
    gaps =  pivot_df.values 
    lons = pivot_df.columns.values 
    lats = pivot_df.index.values 

    datagrid =  xarray.DataArray(data=gaps,
                dims=["lats", "lons"],
                coords={"lons":lons,
                        "lats":lats},
                attrs=dict( description="Gap",
                            units="°",
                        )
                )
    print(datagrid)
    return datagrid,off_stations

def plot_meshgrid(grid,shapes=None,stations=None,
            cmap_args=None,save=None):
    reg = [-73.55, -73.25, 3.95, 4.2]
    fig = pygmt.Figure()
    fig.coast(region=reg,
            projection='M4.5i',
            shorelines=True,
            water='lightblue',
            borders=2,
            rivers="a/blue",
            # land='grey',
            frame=["WNse", "xaf+lx-axis", "yaf+ly-axis"])

    if cmap_args == None:
        cmap_args["path"] = None
    else:
        if os.path.isfile(cmap_args["path"]) == False: 
            pygmt.makecpt(cmap="viridis", 
                            series=cmap_args["limits"],
                            output=cmap_args["path"])

    
    fig.grdimage(grid=grid, cmap=cmap_args["path"])

    # fig.colorbar(frame=["x+lGAP", "y+lm"],  ## barra a la izquierda
    #             position="JML+o1.5c/0c+w8c")
    fig.colorbar(frame=["x+lGAP"], cmap=cmap_args["path"])

    if shapes != None:
        for shp in shapes:
            fig.plot(**shp.__dict__)
            # fig.plot(data=shp.data,color=shp.color)

    if stations != None:
        for status,(sta,lon,lat) in stations.items():
            if sta.size != 0: 
                if status.upper() == "OFF":
                    fig.plot(lon,lat,
                        style='t0.15i',
                        color='red',
                        label='Estaciones',
                        )
                elif status.upper() == "ON":
                    fig.plot(lon,lat,
                            style='t0.15i',
                            color='black',
                            label='Estaciones',
                            )
                xname_padd = abs(reg[0]-reg[1])/40
                yname_padd = abs(reg[2]-reg[3])/40
                fig.text(textfiles=None,x=lon-xname_padd, y=lat+yname_padd, position=None,text=sta,
                            angle=0, font='8p,Helvetica-Bold,black', justify='LM')
    print(fig)
    with pygmt.config(FONT_TITLE=5):
        # fig.basemap(rose="jTL+w1.3c+lO,E,S,N+o-0.1c/3c", map_scale="jBL+w50k+o0.5c/0.5c+f")
        fig.basemap(map_scale="jBR+w10k+o0.5c/0.5c+f+lkm")
    if save == None:
        fig.show()
    else:
        # plt.show()
        if os.path.isdir(os.path.dirname(save)) == False:
            os.makedirs(os.path.dirname(save))
        fig.savefig(save)

def plot_map(grid,shapes=None,stations=None,save=None):
    fig = pygmt.Figure()
    
    # reg1 = [-83, -68, -5, 15] Colombia 
    reg1 = [-85, -71, 2, 13]
    reg2 = [-73.55, -73.25, 3.95, 4.2]

    # reg2 = [-74, -73, 3.5, 4.5]


    # fig.grdimage('@earth_relief_03s',
    fig.grdimage('@earth_relief_05m',
                region=reg1,
                projection='M8i',
                cmap='etopo1',
                shading=True)
    fig.coast(region=reg1,
            projection='M8i',
            shorelines=True,
            borders=['2/1p,black,-','1/2.1p,black'], 
            water='lightblue',
            rivers="a/blue",
            # land='grey',
            
            frame='f')

    fig.plot([reg2[0], reg2[1], reg2[1], reg2[0], reg2[0]],
                [reg2[2], reg2[2], reg2[3], reg2[3], reg2[2]],
        pen="2p,red,-"
        )

    fig.shift_origin(xshift='0.01i',yshift='h-3.75i')  # Shift for next call
    

    # # Second figure
    fig.coast(region=reg2,
            projection='M4.5i',
            shorelines=True,
            water='lightblue',
            borders=2,
            rivers="a/blue",
            # land='grey',
            frame=["WNse", "xaf+lx-axis", "yaf+ly-axis"])
    fig.grdimage(grid=grid, cmap="viridis")

    # fig.colorbar(frame=["x+lGAP", "y+lm"],  ## barra a la izquierda
    #             position="JML+o1.5c/0c+w8c")
    fig.colorbar(frame=["x+lGAP"])

    if shapes != None:
        for data_shp in shapes:
            print(data_shp)
            fig.plot(data=data_shp,color="green")

    if stations != None:
        sta,lon,lat = stations
        fig.plot(lon,lat,
                style='t0.15i',
                color='black',
                label='Estaciones',
                )

    if save == None:
        fig.show()
    else:
        if os.path.isdir(os.path.dirname(save)) == False:
            os.makedirs(os.path.dirname(save))
        fig.savefig(save)

def run_itermap(gap_folder,stations_path,
                shapes=None,cmap_args=None):
    for dp, dn, filenames in os.walk(gap_folder):
        for f in filenames:

            gap_path = os.path.join(dp, f)
            _, file_extension = os.path.splitext(gap_path)
            if file_extension == ".csv":
                print(gap_path)
                pre, ext = os.path.splitext(gap_path)
                save = os.path.join(gap_path, pre + '.png')

                grid,off_stations = get_meshgrid(gap_path)
                print(off_stations )
                stations = get_stations(stations_path,off_stations)
                plot_meshgrid(grid,shapes,stations,cmap_args,save)
                print(save,"\n")
                # print(filepath)

if __name__ == "__main__":
    print("plot_isogap")
    ################### plot meshgrid

    # gap_path = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/outs/6_stations/net_1_S0.005.csv" 
    # stations_path = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/data/stations.csv"
    # save = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/outs/6_stations/net_1_S0.005.png"
    # stations = get_stations(stations_path)
    # grid = get_meshgrid(gap_path)
    # plot_meshgrid(grid,None,stations,save)


    ###################  ITER PLOT
    bloque_shp = "/mnt/SharedDrives/Ecopetrol/qgis/Campos/Apiay/gap/apiay_files/bloque_apiay.shp"
    sh = get_shape(bloque_shp,{"AREA_NOMBR":["APIAY"]})
    sh_obj = ShapeObject(data=sh[0],label="inyeccion",pen="1p,red,.")
    print(sh_obj.__dict__)

    gap_folder = "/home/emmanuel/Ecopetrol/ISOGAP/outs"
    stations_path = "/home/emmanuel/Ecopetrol/ISOGAP/data/stations.csv"
    cmap_args = {"path": "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/data/viridis_80_360.cpt",
                "limits":[80,360]}
    run_itermap(gap_folder,stations_path,shapes=[sh_obj ],cmap_args=cmap_args)
    # run_itermap(gap_folder,stations_path,shapes=None,cmap_args=cmap_args)

    ###################  SHP PROVES
    # bloque_shp = "/mnt/SharedDrives/Ecopetrol/qgis/data/MAPA/MAPA/head shape/Heads-Trajectory.shp"
    # df = gpd.read_file(bloque_shp)
    # sh = get_shape(bloque_shp,{"Name":df["Name"].to_list()})
    # print(sh)
    
    # fig = pygmt.Figure()
    # fig.basemap(region=[-73.55, -73.25, 3.95, 4.2], projection="M4i", frame=True)
    # fig.coast( 
    #     water='skyblue', 
    #     shorelines=True)
    # for  shape in sh:
    #     if len(shape) == 1:
    #         sh_obj = ShapeObject(data=shape,label="inyeccion",color="black")
    #     else:
    #         sh_obj = ShapeObject(data=shape,label="inyeccion",pen="1p,black")
    #     fig.plot(**sh_obj.__dict__)
    # fig.savefig('map1.png')


    # import matplotlib.pyplot as plt

    # plt.plot(range(10))

    # plt.title('Center Title')
    # plt.title('Left Title', loc='left')
    # plt.title('Right Title', loc='right')

    # plt.savefig("/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/prove.png")