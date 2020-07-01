__all__ = [
    "read_IGRF13_COF",
    "read_IGRF13coeffs",
    "read_WMM",
    "read_fortran_DATA",
    "read_gauss_coeff",
    "read_WWW_test_2020",
]


def read_gauss_coeff(file=None):

    if file is None:
        file = "IGRF13.COF"

    if file == "IGRF13.COF":
        (
            dic_dic_h,
            dic_dic_g,
            dic_dic_SV_h,
            dic_dic_SV_g,
            dic_N,
            Years,
        ) = read_IGRF13_COF(file)
    elif file == "IGRF13coeffs.txt":
        (
            dic_dic_h,
            dic_dic_g,
            dic_dic_SV_h,
            dic_dic_SV_g,
            dic_N,
            Years,
        ) = read_IGRF13coeffs(file)
    elif file == "WMM_2015.COF":
        dic_dic_h, dic_dic_g, dic_dic_SV_h, dic_dic_SV_g, dic_N, Years = read_WMM(file)
    elif file == "WMM_2020.COF":
        dic_dic_h, dic_dic_g, dic_dic_SV_h, dic_dic_SV_g, dic_N, Years = read_WMM(file)
    elif file == "FORTRAN_1900_1995.txt":
        dic_dic_h, dic_dic_g, dic_N, Years = read_fortran_DATA(file)
    else:
        raise Exception(f"undefinited file :{file}")

    return dic_dic_h, dic_dic_g, dic_dic_SV_h, dic_dic_SV_g, dic_N, Years


def read_IGRF13_COF(file):

    """read_hg assigns the IGRF13.COF coefficients h and g, in unit of nT, from the text file 
    available along with the Geomag 7.0 software (Windows version) https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
    
    Arguments
        file containing the coefficients h ,g, SVh, SVg versus m,n, year (IGRF13.COF)
    
    Returns
        dic_dic_h (dict of dict): h coefficients {year: {(m,n):h,...},...} year ia string
        dic_dic_g (dict of dict): g coefficients {year: {(m,n):g,...},...} year ia string
        dic_dic_SV_h (dict of dict): SV_h coefficients {year: {(m,n):SV_h,...},...} year ia string
        dic_dic_SV_g (dict of dict): SV_g coefficients {year: {(m,n):SV_g,...},...} year ia string
        dic_N (dict): dictionary containing the order N of the SH decomposition, dic_N[year]=N
        Years (nparray): array of the tabulated year """

    # Standard Library dependencies
    import re
    import os

    # 3rd party dependencies
    import pandas as pd
    import numpy as np

    file = os.path.join(os.path.dirname(__file__), file)

    df = pd.read_table(file, delim_whitespace=True, names=[str(i) for i in range(12)])[
        [str(i) for i in range(6)]
    ]

    indexes_year = [
        (i, x) for i, x in enumerate(list(df["0"])) if ("IGRF" in x) or ("DGRF" in x)
    ]
    indexes = [x[0] for x in indexes_year]
    year = [re.findall(r"\d+", x[1])[0] for x in indexes_year]
    dic_dic_g = {}
    dic_dic_h = {}
    dic_dic_SV_g = {}
    dic_dic_SV_h = {}
    dic_N = {}
    years = []
    dfs = df

    for i, nitems in enumerate(
        np.append(np.diff(indexes), [len(df["0"]) - indexes[-1]])
    ):
        if len(year[i]) == 2:
            year[i] = "19" + year[i]
        dfs = np.split(dfs, [nitems], axis=0)
        dg = dfs[0].iloc[1:]
        dic_dic_g[year[i]] = {
            (int(x[0]), int(x[1])): x[2] for x in zip(dg["1"], dg["0"], dg["2"])
        }
        dic_dic_h[year[i]] = {
            (int(x[0]), int(x[1])): x[2] for x in zip(dg["1"], dg["0"], dg["3"])
        }
        dic_dic_SV_g[year[i]] = {
            (int(x[0]), int(x[1])): x[2] for x in zip(dg["1"], dg["0"], dg["4"])
        }
        dic_dic_SV_h[year[i]] = {
            (int(x[0]), int(x[1])): x[2] for x in zip(dg["1"], dg["0"], dg["5"])
        }
        dic_N[year[i]] = max([x[0] for x in dic_dic_g[year[i]].keys()])
        years.append(float(year[i]))
        dfs = dfs[1]

    years = np.array(years)

    return dic_dic_h, dic_dic_g, dic_dic_SV_h, dic_dic_SV_g, dic_N, np.array(years)


def read_IGRF13coeffs(file):

    """read_hg assigns the IGRF-13 coefficients h and g, in unit of nT, from the text file 
    downloaded from https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
    
    
    Arguments:
        file (string): name of the file (WMM_2015.COF or WMM_2015.COF)
    
    Returns:
        dic_dic_h (dict of dict): h coefficients {year: {(m,n):h,...},...}
        dic_dic_g (dict of dict): g coefficients {year: {(m,n):g,...},...}
        dic_dic_SV_h (dict of dict): SV_h coefficients {year: {(m,n):SV_h,...},...}
        dic_dic_SV_g (dict of dict): SV_g coefficients {year: {(m,n):SV_g,...},...}
        dic_N (dict): dictionary containing the order N of the SH decomposition, dic_N[year]=N
        Years (list): list of the tabulated year """

    # Standard Library dependencies
    import os

    # 3rd party dependencies
    import numpy as np
    import pandas as pd

    file = os.path.join(os.path.dirname(__file__), file)

    df = pd.read_csv(file, header=3, sep="\s+")
    Years = [x for x in df.columns if x[-2:] == ".0"]
    v = []
    for x in df.groupby("g/h"):
        v.append(x)
    g = v[0][1]
    h = v[1][1]
    dic_dic_g = {}
    dic_dic_h = {}
    dic_dic_SV_g = {}
    dic_dic_SV_h = {}
    dic_N = {}
    for Year in Years:
        key_Year = str(int(float(Year)))
        dic_dic_g[key_Year] = {(x[0], x[1]): x[2] for x in zip(g["m"], g["n"], g[Year])}
        dic_dic_h[key_Year] = {(x[0], x[1]): x[2] for x in zip(h["m"], h["n"], h[Year])}
        dic_dic_SV_g[key_Year] = {(x[0], x[1]): 0 for x in zip(g["m"], g["n"])}
        dic_dic_SV_h[key_Year] = {(x[0], x[1]): 0 for x in zip(g["m"], g["n"])}
        index = set([x[0] for x in dic_dic_h[key_Year].keys()])
        N = max(index)
        dic_N[key_Year] = N  # must be 13
        for n in range(1, N + 1):
            dic_dic_h[key_Year][(0, n)] = 0
            dic_dic_SV_h[key_Year][(0, n)] = 0
    dic_dic_SV_h["2020"] = {
        (x[0], x[1]): x[2] for x in zip(h["m"], h["n"], h["2020-25"])
    }
    dic_dic_SV_g["2020"] = {
        (x[0], x[1]): x[2] for x in zip(g["m"], g["n"], g["2020-25"])
    }
    Years = np.array([float(x) for x in Years])
    return dic_dic_h, dic_dic_g, dic_dic_SV_h, dic_dic_SV_g, dic_N, Years


def read_WMM(file):

    """read_hg assigns the WMM coefficients h and g, in unit of nT, from the text file 
    downloaded from https://www.ngdc.noaa.gov/geomag/WMM/wmm_ddownload.shtml
    
    Arguments:
        file (string): name of the file (WMM_2015.COF or WMM_2015.COF)
    
    Returns:
        dic_dic_h (dict of dict): h coefficients {year: {(m,n):h,...},...}
        dic_dic_g (dict of dict): g coefficients {year: {(m,n):g,...},...}
        dic_dic_SV_h (dict of dict): SV_h coefficients {year: {(m,n):SV_h,...},...}
        dic_dic_SV_g (dict of dict): SV_g coefficients {year: {(m,n):SV_g,...},...}
        dic_N (dict): dictionary containing the order N of the SH decomposition, dic_N[year]=N
        Years (list): list of the year povided
    """

    # Standard Library dependencies
    import re
    import os

    # 3rd party dependencies
    import pandas as pd
    import numpy as np

    file = os.path.join(os.path.dirname(__file__), file)

    df = pd.read_csv(file, sep="\s+", skipfooter=2, engine="python")
    df = df.reset_index(level=[0, 1])
    df = df.reset_index()
    df.columns = ["g", "n", "m", "h", "SVg", "SVh"]

    year = re.findall(r"\d+", os.path.basename(file))[0]
    dic_dic_h = {year: {(x[0], x[1]): x[2] for x in zip(df["m"], df["n"], df["h"])}}
    dic_dic_g = {year: {(x[0], x[1]): x[2] for x in zip(df["m"], df["n"], df["g"])}}
    dic_dic_SV_h = {
        year: {(x[0], x[1]): x[2] for x in zip(df["m"], df["n"], df["SVh"])}
    }
    dic_dic_SV_g = {
        year: {(x[0], x[1]): x[2] for x in zip(df["m"], df["n"], df["SVg"])}
    }
    dic_N = {year: max(set([x[0] for x in dic_dic_h[year].keys()]))}
    Years = np.array([float(year)])
    return dic_dic_h, dic_dic_g, dic_dic_SV_h, dic_dic_SV_g, dic_N, Years


def read_fortran_DATA(file):

    """read_hg assigns the coefficients h and g, in unit of nT as 
    extracteed from the FORTRAN program IGRF13 https://www.ngdc.noaa.gov/IAGA/vmod/igrf13.f
    
    Arguments:
        file (string): name of the file (FORTRAN_1900_1995.txt or FORTRAN_2000_2020.txt)
    
    Returns:
        dic_dic_h (dict of dict): h coefficients {year: {(m,n):h,...},...}
        dic_dic_g (dict of dict): g coefficients {year: {(m,n):g,...},...}
        dic_dic_SV_h (dict of dict): SV_h coefficients {year: {(m,n):SV_h,...},...}
        dic_dic_SV_g (dict of dict): SV_g coefficients {year: {(m,n):SV_g,...},...}
        dic_N (dict): dictionary containing the order N of the SH decomposition, dic_N[year]=N
        Years (list): list of the year povided
    """

    # Standard Library dependencies
    import re
    import os

    # 3rd party dependencies
    import pandas as pd
    import numpy as np

    file = os.path.join(os.path.dirname(__file__), file)

    def construct_dic(df):
        df[0] = df[0].apply(lambda x: float(x.split(" ")[-1]))
        df = df.drop([df.columns[-1]], axis=1)
        df = df.T
        res = []
        for x in df.columns:
            res = res + list(df[x])
        N = 0
        while len(res) - N * N - 2 * N > 0:
            N += 1
        N -= 1
        dic_g = {}
        dic_h = {}
        idx = 0
        for n in range(1, N + 1):
            for m in range(0, n + 1):
                dic_g[(m, n)] = res[idx]
                idx += 1
                if m == 0:
                    dic_h[(0, n)] = 0
                else:
                    dic_h[(m, n)] = res[idx]
                    idx += 1
        return N, dic_h, dic_g

    df = pd.read_csv(file, sep=",", header=None, skipfooter=0, engine="python")
    dic_dic_h = {}
    dic_dic_g = {}
    dic_N = {}
    Years = []
    for dg in df.groupby(df.columns[-1]):
        year = str(dg[0])
        dh = dg[1].copy()
        N, dic_h, dic_g = construct_dic(dh)

        dic_dic_h[year] = dic_h
        dic_dic_g[year] = dic_g
        dic_N[year] = N
        Years.append(float(year))
    Years = np.array(Years)

    return dic_dic_h, dic_dic_g, dic_N, Years
    
def read_WWW_test_2020(index):

    """reads  the Test Values for WMM2020 .xlsx file
    
    Arguments:
        index (int): index>=0 and index < 11
    
    Returns:
        Date (dict): Date 
        height (float): height in meters
        colatitude (float): colatitude in ° 
        longitude (float):longitude in ° 
        WMM (dict): 
    """
    
    import pandas as pd
    import os
    
    file = os.path.join(os.path.dirname(__file__), 'WMM2020testvalues.xlsx')

    df = pd.read_excel(file, header=1)
    WMM = df.to_dict()
    WMM = {key:value[index] for key, value in WMM.items()} 
    Date = {"mode":"dec","year":WMM['Date'] }
    height = WMM['Height\n(km)']*1000
    colatitude = 90 - WMM['Lat\n(Deg)']
    longitude = WMM['Lon\n(Deg)']

    del WMM['Date']
    del WMM['Height\n(km)']
    del WMM['Lat\n(Deg)']
    del WMM['Lon\n(Deg)']
    
    return Date, height, colatitude, longitude, WMM
