# /**
#  * @author Emmanuel Castillo
#  * @email ecastillot@unal.edu.co
#  * @create date 2021-08-19 09:33:05
#  * @modify date 2021-08-19 09:33:05
#  * References:
        # https://stackoverflow.com/questions/36556890/restructuring-pandas-dataframe-into-meshgrid-for-basemap (basemap)
        # https://www.pygmt.org/dev/gallery/images/grdclip.html#sphx-glr-gallery-images-grdclip-py (pygmt)
#  */
import os
import pygmt
import xarray
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import geopandas as gpd

def get_shape(shp_path,name_list=None):
    gdf = gpd.read_file(shp_path)
    print(gdf)

    if name_list != None:
        all_data = []
        for name in name_list:
            all_data.append( gdf[gdf["Name"]==name] )
    else:
        all_data = [gdf]
    return all_data

def get_stations(csv_path):
    df = pd.read_csv(csv_path)
    sta = df["STA"].to_numpy()
    lon = df["LON"].to_numpy()
    lat = df["LAT"].to_numpy()
    return sta,lon,lat 

def get_meshgrid(csv_path):

    df = pd.read_csv(csv_path)
    print(df)
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
    return datagrid

def plot_meshgrid(grid,shapes=None,stations=None,save=None):
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
    print(fig)

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

def run_itermap(gap_folder,stations_path):
    for dp, dn, filenames in os.walk(gap_folder):
        for f in filenames:
            gap_path = os.path.join(dp, f)

            pre, ext = os.path.splitext(gap_path)
            save = os.path.join(gap_path, pre + '.png')

            grid = get_meshgrid(gap_path)
            stations = get_stations(stations_path)
            plot_meshgrid(grid,None,stations,save)
            print(save)
            # print(filepath)


if __name__ == "__main__":
    # gap_path = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/outs/6_stations/net_1_S0.005.csv" 
    # stations_path = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/data/stations.csv"
    # save = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/outs/6_stations/net_1_S0.005.png"
    # stations = get_stations(stations_path)
    # grid = get_meshgrid(gap_path)
    # plot_meshgrid(grid,None,stations,save)

    gap_folder = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/outs"
    # gap_folder = "/home/emmanuel/5_stations"
    stations_path = "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/data/stations.csv"
    run_itermap(gap_folder,stations_path)


    # import matplotlib.pyplot as plt

    # plt.plot(range(10))

    # plt.title('Center Title')
    # plt.title('Left Title', loc='left')
    # plt.title('Right Title', loc='right')

    # plt.savefig("/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/prove.png")