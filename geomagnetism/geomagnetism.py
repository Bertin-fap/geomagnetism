__all__ = [
    "B_components",
    "geodetic_to_geocentric",
    "geodetic_to_geocentric_IGRF13",
    "decdeg2dms",
    "Norm_Schimdt",
    "Norm_Stacey",
    "grid_geomagnetic",
    "construct_xarray",
]

# Ellipsoid parameters: semi major axis in metres, reciprocal flattening.
GRS80 = (6_378_137, 298.257222100882711)
WGS84 = (6_378_137, 298.257223563)

# Ellipsoid parameters: semi major axis in metres, semi minor axis in metres.
GRS80_ = (6378.137, 6356.752314140355847852106)
WGS84_ = (6378.137, 6356.752314245179497563967)

# mean radius of the earth in metres must be chosen in accorance with the radius used in determining
# the SH coefficients
EARTH_RADIUS = 6_371_200


def B_components(
    phi_,
    theta_,
    altitude,
    Date,
    referential="geodetic",
    file="IGRF13.COF",
    ELLIPSOID=WGS84,
    SV=False,
):

    """B_components computes the geomagnetic magnetic field components.
    We use the Peddie's notation Peddie, N. W. (1982). "International Geomagnetic Reference Field : the third generation."
    J. Geomag. Geolectr 34: 309-326.
    
    Arguments
        phi_ (float): longitude in deg (0° - 360°) or West East longitude (-180° - +180°)
        theta_(float): Colatitude in deg (0° - 180°) 
        altitude (float): Elevation in m,
        Date :  Date used to compute the magnetic field 1900<= Date< 2025 . Date is a dictionary
                Date["mode"]="ymd" if Date is expressed as yyyy/mm/dd otherwise Date is in decimal value
                Date["year"]= year if Date["mode"]="ymd" oterwise Date["year"]=decimal value year
                Date["month"]
                Date["day"],
                Date["hour"]
                Date["minute"]
                Date["second"]
            ex: Date = {"mode":"ymd","year":2020,"month":16,"day":16,"hour":0,"minute":0,"second":0}
                Date = {"dec":"ymd","year":2020.5}
        referential : referential = geotetic  if the colatitude is expressed in geotetic referential
                    : referential = geocentric if the colatitude is expressed in geotetic referential
        file (string) : name of the file containing the Gauss coefficients h and g 
                      file = "IGRF13.COF" or "WMM_2015.COF" or "WMM_2020.COF" (default "WMM_2020.COF")
        ELLIPSOID (tuple): WGS84 or WGS84 (default WGS84)
        SV (bool): if True computation of the field secular variation
        
    Returns
    result (dict):{'D': Declination in deg,
                   'F': Total field intensity in nT,
                   'H': Horizontal field intensity in nT,
                   'I': Inclination in deg,
                   'X': North component in nT in the geocentric coordinate,
                   'Y': East component in nT in both geocentric and geotetic coordinate,
                   'Z': Down component in nT in the geocentric coordinate,
                   'Fd': Total field intensity secular variation in nT/year if SV==True None elsewhere,
                   'Hd': Horizontal field intensity secular variation in nT/year if SV==True None  if SV==True None elsewhere,
                   'Yd': East component secular variation in nT/year in both geocentric and geotetic coordinate if SV==True None elsewhere,
                   'Zd': Down component secular variatin in nT/year in the geocentric coordinate if SV==True None elsewhere,
                   'Id': Inclination secular variationin deg/year,
                   'Dd': Declination secular variationin deg/year
    """

    # Standard Library dependencies
    import os as os

    # 3rd party import
    import numpy as np
    from scipy.special import lpmn

    # Internal import
    from .read_geo_data import read_gauss_coeff

    EPS = 1.0e-5  # if colatitude <EOS rad or colatitude > pi-EPS we call field_computation_pole
    ADD_YEARS = 5  # extrapolation possible over ADD_YEARS additional years

    def assign_hg(Date):

        """Time interpollation/extrapollation of the coefficients h and g :
        - for years[0]<= Year <years[-1] we use an interpollation scheme
        - for years[-1]<=year<=years[-1]+ 5 we use an extrapolationllation scheme with a secular variation (SV) nT/year """

        if Date["mode"] == "ymd":
            del Date["mode"]
            Year = decimal_year(**Date)
        else:
            Year = Date["year"]

        if (Year < years[0]) or (Year > years[-1] + ADD_YEARS):
            raise Exception(
                f"Year={Year} out of range. Must be >={years[0]} and <= {years[0]+5}"
            )

        idx = (np.where(years <= Year))[0][-1]
        year = int(years[idx])
        dt = Year - year
        dic_h0 = dic_dic_h[str(year)]
        dic_g0 = dic_dic_g[str(year)]
        N = dic_N[str(year)]
        dic_h = {}
        dic_g = {}
        if (
            (Year >= years[0]) & (Year < years[-1]) & (years[0] != years[-1])
        ):  # use interpollation
            Dt = years[idx + 1] - years[idx]
            dic_h1 = dic_dic_h[str(int(years[idx + 1]))]
            dic_g1 = dic_dic_g[str(int(years[idx + 1]))]
            for x in dic_h0.keys():
                dic_h[x] = dic_h0[x] + dt * (dic_h1[x] - dic_h0[x]) / Dt
                dic_g[x] = dic_g0[x] + dt * (dic_g1[x] - dic_g0[x]) / Dt

        elif (Year >= years[-1]) & (
            Year < years[-1] + ADD_YEARS
        ):  # use extrapolation through secular variation (SV)
            for x in dic_h0.keys():
                dic_h[x] = dic_h0[x] + dic_dic_SV_h[str(int(years[-1]))][x] * dt
                dic_g[x] = dic_g0[x] + dic_dic_SV_g[str(int(years[-1]))][x] * dt

        return dic_h, dic_g, N, year

    # Assign dic_dic_g,dic_dic_h,dic_dic_SV_g,dic_dic_SV_h,dic_N,years
    # Avoid the reaffectation of these coefficients after the first call of B_components

    B_components.flag = getattr(B_components, "flag", True)
    if B_components.flag:  # First call
        (
            dic_dic_h,
            dic_dic_g,
            dic_dic_SV_h,
            dic_dic_SV_g,
            dic_N,
            years,
        ) = read_gauss_coeff(file=file)
        B_components.dic_dic_g = getattr(B_components, "dic_dic_g", dic_dic_g)
        B_components.dic_dic_h = getattr(B_components, "dic_dic_h", dic_dic_h)
        B_components.dic_dic_SV_h = getattr(B_components, "dic_dic_SV_h", dic_dic_SV_h)
        B_components.dic_dic_SV_g = getattr(B_components, "dic_dic_SV_g", dic_dic_SV_g)
        B_components.dic_N = getattr(B_components, "dic_N", dic_N)
        B_components.years = getattr(B_components, "years", years)
        B_components.flag = False

    else:
        dic_dic_g = B_components.dic_dic_g
        dic_dic_h = B_components.dic_dic_h
        dic_dic_SV_h = B_components.dic_dic_SV_h
        dic_dic_SV_g = B_components.dic_dic_SV_g
        dic_N = B_components.dic_N
        years = B_components.years

    dic_h, dic_g, N, year = assign_hg(
        Date
    )  # performs interpolation of the coefficients h and g

    # Compute the transformation matrix from geodetic to geocentric frames
    if referential == "geodetic":
        r_geocentric, co_latitude_geocentric, delta = geodetic_to_geocentric(
            ELLIPSOID, theta_, altitude
        )
        theta = co_latitude_geocentric
        mat_rot = np.array(
            [
                [np.cos(delta), 0, np.sin(delta)],
                [0, 1, 0],
                [-np.sin(delta), 0, np.cos(delta)],
            ]
        )
    else:
        r_geocentric = EARTH_RADIUS + altitude
        theta = theta_ * np.pi / 180
        mat_rot = np.identity(3)

    r_a = EARTH_RADIUS / r_geocentric

    if phi_ >= 0:
        phi = phi_ * np.pi / 180
    else:
        phi = (360 + phi_) * np.pi / 180

    # synthesis of Xc, Yc and Zc in geocentric coordinates
    if theta > EPS and theta < np.pi - EPS:
        # Legendre polynomes computation
        Norm = Norm_Schimdt(N, N)
        M, Mp = lpmn(N, N, np.cos(theta))
        M = M * Norm
        Mp = Mp * Norm * (-1) * np.sin(theta)  # dPn,m(cos(theta))/d theta
        X, Y, Z = field_computation(r_a, M, Mp, phi, theta, dic_g, dic_h, N, mat_rot)

    else:
        X, Y, Z = field_computation_pole(r_a, phi, theta, dic_g, dic_h, N, mat_rot, EPS)

    F = np.linalg.norm([X, Y, Z])
    H = np.linalg.norm([X, Y])
    D = np.arctan2(Y, X) * 180 / np.pi  # declination
    I = np.arctan2(Z, H) * 180 / np.pi  # inclination

    # secular variation computation
    Xd, Yd, Zd, Hd, Fd = [None] * 5
    if SV:
        if theta > EPS and theta < np.pi - EPS:
            # Legendre polynomes computation
            Norm = Norm_Schimdt(N, N)
            M, Mp = lpmn(N, N, np.cos(theta))
            M = M * Norm
            Mp = Mp * Norm * (-1) * np.sin(theta)  # dPn,m(cos(theta))/d theta
            Xd, Yd, Zd = field_computation(
                r_a,
                M,
                Mp,
                phi,
                theta,
                dic_dic_SV_g[str(year)],
                dic_dic_SV_h[str(year)],
                N,
                mat_rot,
            )

        else:
            Xd, Yd, Zd = field_computation_pole(
                r_a,
                phi,
                theta,
                dic_dic_SV_g[str(year)],
                dic_dic_SV_h[str(year)],
                N,
                mat_rot,
                EPS,
            )

        Hd = (X * Xd + Y * Yd) / H
        Fd = (X * Xd + Y * Yd + Z * Zd) / F
        Id = 180 * (H * Zd - Z * Hd) / (F * F * np.pi)
        Dd = 180 * (X * Yd - Y * Xd) / (H * H * np.pi)

    result = dict(
        zip(
            [
                "X",
                "Y",
                "Z",
                "F",
                "H",
                "I",
                "D",
                "Xd",
                "Yd",
                "Zd",
                "Hd",
                "Fd",
                "Id",
                "Dd",
            ],
            [X, Y, Z, F, H, I, D, Xd, Yd, Zd, Hd, Fd, Id, Dd],
        )
    )

    return result


def field_computation_pole(r_a, phi, theta, dic_g, dic_h, N, mat_rot, EPS):

    """
    compute the geomagnetic field in the geotetic coordinates near the north and south pole.
    
    Arguments:
        r_a (float): geocentric radial ccordinale/radius of the Earth
        phi (float); longitude (rad)
        theta: 
        dic_g (dict): g Gauss coefficients 
        dic_h (dict): h Gauss coefficients
        N (int): order of the spherical harmonic decomposition
    Returns:
        X (float): North component in nT in the geocentric coordinate,
        Y (float): East component in nT in both geocentric and geotetic coordinate,
        Z (float): Down component in nT in the geocentric coordinate,
    
    """
    import numpy as np

    Zc = 0
    Sh = 0
    Sg = 0
    coef = r_a * r_a
    if theta < EPS:
        for n in range(1, N + 1):
            coef *= r_a
            Zc -= (n + 1) * (r_a) ** (n + 2) * dic_g[(0, n)]
            Sh += np.sqrt(n * (n + 1) / 2) * (r_a) ** (n + 2) * dic_h[(1, n)]
            Sg += np.sqrt(n * (n + 1) / 2) * (r_a) ** (n + 2) * dic_g[(1, n)]

        Xc = Sg * np.cos(phi) + Sh * np.sin(phi)
        Yc = Sg * np.sin(phi) - Sh * np.cos(phi)

    elif theta > np.pi - EPS:
        for n in range(1, N + 1):
            coef *= -r_a
            Zc -= coef * (n + 1) * dic_g[(0, n)]
            Sh += coef * np.sqrt(n * (n + 1) / 2) * dic_h[(1, n)]
            Sg += coef * np.sqrt(n * (n + 1) / 2) * dic_g[(1, n)]

        Xc = Sg * np.cos(phi) + Sh * np.sin(phi)
        Yc = -Sg * np.sin(phi) + Sh * np.cos(phi)

    [X, Y, Z] = np.matmul(mat_rot, [Xc, Yc, Zc])

    return X, Y, Z


def field_computation(r_a, M, Mp, phi, theta, dic_g, dic_h, N, mat_rot):
    """
    compute the geomagnetic field in the geotetic coordinates
    Arguments:
        r_a (float): geocentric radial ccordinale/radius of the Earth
        M (np array): matrix of the Legendre polynomes values
        Mp (np array): matrix of the derivative of the Legendre polynomes vs the colatitude
        phi (float); longitude (rad)
        theta: 
        dic_g (dict): g Gauss coefficients 
        dic_h (dict): h Gauss coefficients
        N (int): order of the spherical harmonic decomposition
    Returns:
        X (float): North component in nT in the geocentric coordinate,
        Y (float): East component in nT in both geocentric and geotetic coordinate,
        Z (float): Down component in nT in the geocentric coordinate,
    
    """
    import numpy as np

    # synthesis of Xc, Yc and Zc in geocentric coordinates
    Xc = 0.0
    Yc = 0.0
    Zc = 0.0
    coef = r_a * r_a
    for n in range(1, N + 1):
        coef *= r_a
        Xc += (
            sum(
                [
                    Mp[m, n]
                    * (
                        dic_g[(m, n)] * np.cos(m * phi)
                        + dic_h[(m, n)] * np.sin(m * phi)
                    )
                    for m in range(0, n + 1)
                ]
            )
            * coef
        )
        Yc += (
            sum(
                [
                    m
                    * M[m, n]
                    * (
                        dic_g[(m, n)] * np.sin(m * phi)
                        - dic_h[(m, n)] * np.cos(m * phi)
                    )
                    for m in range(0, n + 1)
                ]
            )
            * coef
        )
        Zc -= (
            sum(
                [
                    M[m, n]
                    * (
                        dic_g[(m, n)] * np.cos(m * phi)
                        + dic_h[(m, n)] * np.sin(m * phi)
                    )
                    for m in range(0, n + 1)
                ]
            )
            * coef
            * (n + 1)
        )
    Yc = Yc / np.sin(theta)

    # conversion to the coordinate system specified by variable referential
    [X, Y, Z] = np.matmul(mat_rot, [Xc, Yc, Zc])

    return X, Y, Z


def geodetic_to_geocentric(ellipsoid, co_latitude, height):

    """Return geocentric (Cartesian) radius and colatitude corresponding to
    the geodetic coordinates given by colatitude (in
    degrees ) and height (in metre) above ellipsoid. 
    see http://clynchg3c.com/Technote/geodesy/coordcvt.pdf
    credit : https://codereview.stackexchange.com/questions/195933/convert-geodetic-coordinates-to-geocentric-cartesian
    with minor modifications.
    
    Arguments:
        ellipsoid (tuple): ellipsoid parameters (semi-major axis, reciprocal flattening)
        co_latitude (float): geotetic colatitude (in degrees) 0°<=co_latitude<=180°
        height (float): height (in metre) above ellipsoid
    
    Returns
        r_geocentric (float): geocentric radius (m)
        co_latitude_geocentric (float): geocentric colatitude (radians)
        delta (float): angle between geocentric and geodetic referentials (radians)
        """

    # 3rd Party dependencies
    from math import radians
    import numpy as np

    lat = radians(90 - co_latitude)  # geodetic latitude
    sin_lat = np.sin(lat)
    a, rf = ellipsoid  # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    r_n = a / np.sqrt(1 - e2 * sin_lat ** 2)  # prime vertical curvature radius
    r = (r_n + height) * np.cos(lat)  # perpendicular distance from z axis
    z = (r_n * (1 - e2) + height) * sin_lat
    r_geocentric = np.sqrt(r ** 2 + z ** 2)
    co_latitude_geocentric = np.pi / 2 - np.arctan(
        (1 - e2 * r_n / (r_n + height)) * np.tan(lat)
    )  # geocentric colatitude
    delta = co_latitude_geocentric - radians(
        co_latitude
    )  # angle between geocentric and geodetic referentials

    return r_geocentric, co_latitude_geocentric, delta


def geodetic_to_geocentric_IGRF13(ellipsoid, co_latitude, height):

    """conversion from geodetic to geocentric coordinates. 
    Translated from FORTRAN program IGRF13
    ellipsoid = GRS80 or WGS84 according to the choice of world geodetic system
    [1] Peddie, Norman W. International Geomagnetic Reference Field : the third generation J. Geomag. Geolectr 34 p. 309-326
    
    Arguments
        ellipsoid (tuple): (a=semi major axis in metres, b=semi major axis in metres)
        colatitude : geodetic colatitude in degree (0<= colatitude <=180)
        height : elevation in meters above the geoid
    Returns
        r : geocentric radius
        ct: cos(geocentric colatitude )
        st: sin(geocentric colatitude )
        cd: cos(delta ) delta is the angle between geocentric and geodetic colatitude (radians)
        sd: sin(delta )"""

    # 3rd Party dependencies
    from math import radians
    import numpy as np

    theta = radians(co_latitude)
    st = np.sin(theta)
    ct = np.cos(theta)
    a, b = ellipsoid  # a,b semi major and semi minor axis
    a2 = a * a
    b2 = b * b
    one = a2 * st * st
    two = b2 * ct * ct
    three = one + two
    rho = np.sqrt(three)  # a*a/r_n with r_n the prime vertical curvature radius
    r = np.sqrt(
        height * (height + 2.0 * rho) + (a2 * one + b2 * two) / three
    )  # geocentric radius [1](6)
    cd = (height + rho) / r  # cos(delta )
    sd = (a2 - b2) / rho * ct * st / r  # sin(delta )
    one = ct
    ct = ct * cd - st * sd  # cos(geocentric colatitude )
    st = st * cd + one * sd  # sin(geocentric colatitude )
    return r, ct, st, cd, sd


def Norm_Schimdt(m, n):

    """Norm_Schimdt buids the normalization matrix which coefficients can be found in
    Winch D. E. et al. Geomagnetism and Schmidt quasi-normalization Geophys. J. Int. 160 p. 487-454 2005"""

    # 3rd party dependencies
    import math
    import numpy as np

    sgn = lambda m: 1 if m % 2 == 0 else -1
    norm = (
        lambda m, n: sgn(m)
        * np.sqrt((2 - (m == 0)) * math.factorial(n - m) / math.factorial(n + m))
        if m > 0
        else 1
    )
    norm_schimdt = []
    for m_ in range(m + 1):
        norm_schimdt.append(
            [norm(m_, n_) if (n_ - np.abs(m_) >= 0) else 0 for n_ in range(n + 1)]
        )
    return np.array(norm_schimdt)


def Norm_Stacey(m, n):

    """Norm_Stacey buids the normalization matrix which coefficients can be found in
    Stracey F. D. et al. Physics of the earth Cambridge appendix C"""

    # 3rd party dependencies
    import math
    import numpy as np

    sgn = lambda m: 1 if m % 2 == 0 else -1
    norm = (
        lambda m, n: sgn(m)
        * np.sqrt(
            (2 - (m == 0)) * (2 * m + 1) * math.factorial(n - m) / math.factorial(n + m)
        )
        if m > 0
        else 1
    )
    norm_stacey = []
    for m_ in range(m + 1):
        norm_stacey.append(
            [norm(m_, n_) if (n_ - np.abs(m_) >= 0) else 0 for n_ in range(n + 1)]
        )
    return np.array(norm_stacey)


def decimal_year(year, month, day, hour, minute, second):

    """decimal_year converts a date (year,month,day,hour,minute,second) into a decimal date. credit : Kimvais
     https://stackoverflow.com/questions/6451655/how-to-convert-python-datetime-dates-to-decimal-float-years"""

    # Standard Library dependencies
    from datetime import datetime

    d = datetime(year, month, day, hour, minute, second)
    year_ = (float(d.strftime("%j")) - 1) / 366 + float(d.strftime("%Y"))

    return year_


def decdeg2dms(dd):

    """ Tansform decimal degrees into degrees minutes seconds
    
    Argument:
        dd (float): decimal angle 
    Returns:
        degrees, minutes, seconds"""

    negative = dd < 0
    dd = abs(dd)
    minutes, seconds = divmod(dd * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    if negative:
        if degrees > 0:
            degrees = -degrees
        elif minutes > 0:
            minutes = -minutes
        else:
            seconds = -seconds
    return degrees, minutes, seconds


def construct_xarray(intensities, angles, longitudes, colatitudes):

    """
    construct the xarray of a hyperspectrum

    Arguments:
        intensities (np array): geomagnetic field components intensity (nT)
        angles (np array): geomagnetic field inclination and declination
        longitudes (np array): longitude (°)
        colatitudes (np array): colatitude (°)

    Returns:
        dintensities (xarray): hyperspectrum containing the geomagnetic field components intensity
        dangles (xarray): hyperspectrum containing the geomagnetic field declination and inclination
    """

    # 3rd party dependencies
    import xarray as xr
    import numpy as np

    miller = lambda x: 1.25 * np.arcsinh(np.tan(4 * x * np.pi / 900)) * 180 / np.pi

    y = miller(90 - colatitudes)

    dintensities = xr.DataArray(
        intensities,
        dims=["typ", "lat", "long"],
        name="Field intensity",
        attrs={"units": "nT",},
        coords={
            "typ": xr.DataArray(["X", "Y", "Z", "H", "F",], name="typ", dims=["typ"]),
            "lat": xr.DataArray(y, name="lat", dims=["lat"], attrs={"units": "°"}),
            "long": xr.DataArray(
                longitudes, name="long", dims=["long"], attrs={"units": "°"}
            ),
        },
    )

    dangles = xr.DataArray(
        angles,
        dims=["typ", "lat", "long"],
        name="Angle",
        attrs={"units": "°",},
        coords={
            "typ": xr.DataArray(["D", "I",], name="typ", dims=["typ"]),
            "lat": xr.DataArray(y, name="lat", dims=["lat"], attrs={"units": "°"}),
            "long": xr.DataArray(
                longitudes, name="long", dims=["long"], attrs={"units": "°"}
            ),
        },
    )
    return dintensities, dangles


def grid_geomagnetic(colatitudes, longitudes, height=0, Date={"mode":"dec","year":2020.0}):

    """computes the geomagnetic field characteristics on a mesh of
    colatitudes, latitudes
    
    Arguments:
        colatitudes (nparray): list of colatitudes
        latitudes (nparray): list of latitudes
        height (float): height (meters)
     
    Returns:
        da (xarray): containing the values of the geomagnetic field components, 
    declination and inclination
     
    """

    # 3rd party dependencies
    import numpy as np
    
    H = []
    X = []
    Y = []
    Z = []
    D = []
    I = []
    F = []
    result = []

    for colatitude in colatitudes:
        for longitude in longitudes:
            result = B_components(
                longitude,
                colatitude,
                height,
                Date,
                referential="geodetic",
                file="WMM_2020.COF",
                SV=True,
            )
            X.append(result["X"])
            Y.append(result["Y"])
            Z.append(result["Z"])
            H.append(result["H"])
            F.append(result["F"])
            D.append(result["D"])
            I.append(result["I"])

    intensities = [
        np.array(X).reshape((len(colatitudes), len(longitudes))),
        np.array(Y).reshape((len(colatitudes), len(longitudes))),
        np.array(Z).reshape((len(colatitudes), len(longitudes))),
        np.array(H).reshape((len(colatitudes), len(longitudes))),
        np.array(F).reshape((len(colatitudes), len(longitudes))),
    ]

    angles = [
        np.array(D).reshape((len(colatitudes), len(longitudes))),
        np.array(I).reshape((len(colatitudes), len(longitudes))),
    ]

    dintensities, dangles = construct_xarray(
        intensities, angles, longitudes, colatitudes
    )

    return dintensities, dangles
