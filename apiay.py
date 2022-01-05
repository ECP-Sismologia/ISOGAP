# /**
#  * @author [Emmanuel Castillo]
#  * @email [excastillot@unal.edu.co]
#  * @create date 2021-08-27 23:06:52
#  * @modify date 2021-08-27 23:06:52
#  * @desc [description]
#  */

import geopandas as gpd
from plot_isogap import ShapeObject, get_shape,run_itermap

def get_heads_trajectory():
    bloque_shp = "/mnt/SharedDrives/Ecopetrol/qgis/data/MAPA/MAPA/head shape/Heads-Trajectory.shp"
    df = gpd.read_file(bloque_shp)
    sh = get_shape(bloque_shp,{"Name":df["Name"].to_list()})
    shapes = []
    for  shape in sh:
        if len(shape) == 1:
            sh_obj = ShapeObject(data=shape,label="inyeccion",color="pink")
        else:
            sh_obj = ShapeObject(data=shape,label="inyeccion",pen="0.5p,pink")
        shapes.append(sh_obj)
    return shapes

def get_fields():
    apiay = "/mnt/SharedDrives/Ecopetrol/qgis/data/MAPA/MAPA/Poligonos Campos/Campo_Apiay.shp"
    lr = "/mnt/SharedDrives/Ecopetrol/qgis/data/MAPA/MAPA/Poligonos Campos/Campo_LR.shp"
    pachaquiaro = "/mnt/SharedDrives/Ecopetrol/qgis/data/MAPA/MAPA/Poligonos Campos/CAmpo_Pachaquiaro.shp"
    suria = "/mnt/SharedDrives/Ecopetrol/qgis/data/MAPA/MAPA/Poligonos Campos/Campo_Suria.shp"

    shapes = []
    for shp_field in [apiay,lr,pachaquiaro,suria]:
        sh = get_shape(shp_field)
        for shape in sh:
            shp = ShapeObject(data=shape,label="campos",pen="1p,red,.")
            shapes.append(shp)
    print(shapes)
    # shapes = [item for sublist in shapes for item in sublist]
    return shapes

def get_bloque():
    bloque_shp = "/mnt/SharedDrives/Ecopetrol/qgis/Campos/Apiay/gap/apiay_files/bloque_apiay.shp"
    sh = get_shape(bloque_shp,{"AREA_NOMBR":["APIAY"]})
    sh_obj = ShapeObject(data=sh[0],label="bloque",pen="1p,black,.")
    print(sh_obj)
    return [sh_obj]

if __name__ == "__main__":
    #cambios
    campos = get_fields()
    bloque = get_bloque()
    # heads_traj = get_heads_trajectory()
    # print([campos,bloque])
    shapes = [item for sublist in [campos,bloque] for item in sublist]
    gap_folder = "/home/emmanuel/Ecopetrol/ISOGAP/outs"
    stations_path = "/home/emmanuel/Ecopetrol/ISOGAP/data/stations.csv"
    cmap_args = {"path": "/mnt/SharedDrives/Ecopetrol/opt/ISOGAP/data/viridis_80_360.cpt",
                "limits":[80,360]}
    run_itermap(gap_folder,stations_path,shapes=shapes,cmap_args=cmap_args)