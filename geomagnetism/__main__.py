def plot_isolines():


    # 3rd party imports
    import numpy as np
    
    # Internal imports
    import geomagnetism as geo

    colatitudes = np.linspace(0, 180, 181)
    longitudes = np.linspace(-180, 179, 360)

    dintensities, dangles, dintensities_sv, dangles_sv = geo.grid_geomagnetic(
        colatitudes, longitudes
    )
    geo.plot_geomagetism(
        dintensities,
        dangles,
        dintensities_sv,
        dangles_sv,
        "Id",
        {"proj": "spstere", "boundinglat": -55, "lon_0": 270},
    )

    geo.plot_geomagetism(
        dintensities,
        dangles,
        dintensities_sv,
        dangles_sv,
        "Id",
        {"proj": "npstere", "boundinglat": 70, "lon_0": 270},
    )

    geo.plot_geomagetism(
        dintensities,
        dangles,
        dintensities_sv,
        dangles_sv,
        "Id",
        {
            "proj": "mill",
            "llcrnrlat": -90,
            "urcrnrlat": 90,
            "llcrnrlon": 0,
            "urcrnrlon": 360,
        },
    )


if __name__ == "__main__":
    plot_isolines()
