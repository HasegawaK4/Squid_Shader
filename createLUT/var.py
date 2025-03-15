import sys
from utils import *

over_write_variable_dict = {
    "project_name": "str",
    "LUT_SIZE": "int",
    "ANGLE_MAX": "int",
    "FILMNUM": "int",
    "N_THIN": "float",
    "N_AIR": "float",
    "N_THICK": "float",
    "LMDMIN": "int",
    "DO_CALCULATIONS": "bool",
    "USE_MULTIPROCESS": "bool",
    "SAVE_NDARRAY": "bool",
    "SAVE_CSV": "bool",
    "LOAD_NDARRAY": "bool",
    "LOAD_CSV": "bool",
}

PROJECT_NAME = "test05"

LMDMIN = 390

LUT_SIZE = 16
ANGLE_MAX = 90
FILMNUM = 32
N_THIN = 1.25
N_AIR = 1.0
N_THICK = 1.53
D_THICK = 500   # アラゴナイトクリスタル層(単位:nm)
D_THIN = 25     # タンパク質(コンキオリン)層(単位:nm)

DO_CALCULATIONS = True
USE_MULTIPROCESS = True
SAVE_NDARRAY = True
SAVE_CSV = True
LOAD_NDARRAY = False
LOAD_CSV = False

if len(sys.argv) > 1:
    # コマンドライン引数から "変数名=値" という形式を探す
    for arg in sys.argv[1:]:
        for var in over_write_variable_dict.keys():
            if arg.startswith(var + "="):
                # 値を更新
                value = judge_value(arg, over_write_variable_dict[var], var)
                exec(f"{var} = {value}", globals())

SAVE_NAME = f"{PROJECT_NAME}_{LUT_SIZE}"

PRINTLOG_PATH = make_filename_by_seq("./printlog", f"{SAVE_NAME}.txt")

CMF = csv_read('lin2012xyz2e_5_7sf.csv', LMDMIN)
LIGHT_COL = csv_read('D65.csv', LMDMIN)

N_LIST = [N_AIR]
D_LIST = []
for i in range(FILMNUM + 1):  #最後の層までいれるため+1した回数ループ
    if i % 2 == 0:
        N_LIST.append(N_THICK)
    else:
        N_LIST.append(N_THIN)

for i in range(FILMNUM):
    if i % 2 == 0:
        D_LIST.append(D_THICK)
    else:
        D_LIST.append(D_THIN)