import numpy as np
import csv
from scipy import integrate
from utils import *
from var import *
import time
import multiprocessing as mp
import datetime


def clc_st_color(angle, ref_angle, oblique_angle):
    if float(angle) + oblique_angle < 90:
        angle = float(angle) + oblique_angle
    else:
        angle = 180 - float(angle) - oblique_angle
    ref_angle = np.abs(float(ref_angle) - oblique_angle)
    # print(angle, ref_angle)
    # print(self.cmf)
    wavelengths = CMF[:, 0]
    xyz = np.zeros(3)
    # I = integrate.simpson(cmf[:, 2], x=wavelengths)
    r_list, t_list, phi_list = init_param(angle, ref_angle, wavelengths)
    # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
    f_col = f_color(angle, r_list, t_list, phi_list)
    xyz[0] = integrate.simpson(f_col[:, 0], x=wavelengths)
    xyz[1] = integrate.simpson(f_col[:, 1], x=wavelengths)
    xyz[2] = integrate.simpson(f_col[:, 2], x=wavelengths)
    kk = integrate.simpson(f2(2), x=wavelengths)
    k = 1 / kk
    xyz *= k
    rgb = xyz_to_widegamut(xyz)
    return rgb

def init_param(angle, ref_angle, wavelengths):
    r_list, t_list = calc_rt(angle)
    phi_list = calc_phi(angle, ref_angle, wavelengths)
    return r_list, t_list, phi_list

def calc_rt(angle):
    r_list = []
    t_list = []
    tmp_angle = np.radians(angle)
    for i in range(FILMNUM + 1):  #境界の数はfilmnum+1
        tmp_ref_arc_angle = np.arcsin(np.sin(tmp_angle) * N_LIST[i] / N_LIST[i + 1])
        cos_tmp_angle = np.cos(tmp_angle)
        cos_tmp_ref_arc_angle = np.cos(tmp_ref_arc_angle)

        r_s = (N_LIST[i] * cos_tmp_angle - N_LIST[i + 1] * cos_tmp_ref_arc_angle) / (
                N_LIST[i] * cos_tmp_angle + N_LIST[i + 1] * cos_tmp_ref_arc_angle)
        r_p = (N_LIST[i + 1] * cos_tmp_angle - N_LIST[i] * cos_tmp_ref_arc_angle) / (
                N_LIST[i + 1] * cos_tmp_angle + N_LIST[i] * cos_tmp_ref_arc_angle)
        r = (r_s + r_p) / 2

        t_s = N_LIST[i] * cos_tmp_angle / (
                N_LIST[i] * cos_tmp_angle + N_LIST[i + 1] * cos_tmp_ref_arc_angle)
        t_p = N_LIST[i] * cos_tmp_angle / (
                N_LIST[i + 1] * cos_tmp_angle + N_LIST[i] * cos_tmp_ref_arc_angle)
        t = t_s + t_p

        r_list.append(r)
        t_list.append(t)
        tmp_angle = tmp_ref_arc_angle
    return r_list, t_list

def calc_phi(angle, ref_angle, wavelength):  #air
    phi_list = []
    tmp_angle = np.radians(angle)
    tmp_ref_arc_ref_angle = np.arcsin(np.sin(np.radians(ref_angle)) * N_LIST[0] / N_LIST[1])
    for i in range(FILMNUM):
        tmp_ref_arc_angle = np.arcsin(np.sin(tmp_angle) * N_LIST[i] / N_LIST[i + 1])
        delta = calc_path_diff(N_LIST[i], N_LIST[i + 1], D_LIST[i], tmp_angle, tmp_ref_arc_angle,
                                    tmp_ref_arc_ref_angle)
        phi = 2 * np.pi * delta / wavelength
        phi_list.append(phi)
        tmp_angle = tmp_ref_arc_angle
        tmp_ref_arc_ref_angle = np.arcsin(np.sin(tmp_ref_arc_ref_angle) * N_LIST[i + 1] / N_LIST[i + 2])
    return phi_list

def calc_path_diff(n1, n2, d, angle, ref_arc_angle, ref_arc_ref_angle):  #角度は弧度法
    return n2 * d * (1 / np.cos(ref_arc_angle) + 1 / np.cos(ref_arc_ref_angle)) - n1 * d * (
            np.tan(ref_arc_angle) + np.tan(ref_arc_ref_angle)) * np.sin(angle)  #異角度光路差

def f_color(angle, r_list, t_list, phi_list):
    r = clc_final_reflectance(angle, r_list, t_list, phi_list)
    col = r.reshape(-1, 1) * CMF[:, 1:] * LIGHT_COL[:, 1].reshape(-1, 1)
    return col

def f2(index):
    col = CMF[:, index] * LIGHT_COL[:, 1]
    return col

def clc_final_reflectance(angle, r_list, t_list, phi_list):  #層が複数になった最終的な反射率を求める
    if angle < 90:
        Rlist = [r_list[FILMNUM]]  #gamma_1から順に格納
        #print(self.filmnum)
        for i in range(FILMNUM, 0, -1):
            gamma = (r_list[i - 1] + Rlist[FILMNUM - i] * (
                    np.cos(2 * phi_list[i - 1]) + 1j * np.sin(2 * phi_list[i - 1]))) / (
                            1 + r_list[i - 1] * Rlist[FILMNUM - i] * (
                            np.cos(2 * phi_list[i - 1]) + 1j * np.sin(2 * phi_list[i - 1])))
            Rlist.append(gamma)
        return abs(Rlist[FILMNUM]) ** 2
    else:
        #print(self.t)
        Tlist = [t_list[FILMNUM]]
        for i in range(FILMNUM, 0, -1):
            gamma = (t_list[i - 1] * Tlist[FILMNUM - i] * (
                    np.cos(phi_list[i - 1]) + 1j * np.sin(phi_list[i - 1]))) / (
                            1 + r_list[i - 1] * Tlist[FILMNUM - i] * (
                            np.cos(2 * phi_list[i - 1]) + 1j * np.sin(2 * phi_list[i - 1])))
            Tlist.append(gamma)
        return abs(Tlist[FILMNUM]) ** 2



# 1つのoblique_angleに対する処理を並列化
def process_oblique_angle(oblique_angle):
    result = np.zeros((LUT_SIZE, LUT_SIZE, 3))
    for angle in range(LUT_SIZE):
        print(oblique_angle, angle)
        for ref_angle in range(LUT_SIZE):
            angle_ = angle / (LUT_SIZE - 1) * ANGLE_MAX
            oblique_angle_ = oblique_angle / (LUT_SIZE - 1) * ANGLE_MAX
            ref_angle_ = ref_angle / (LUT_SIZE - 1) * ANGLE_MAX
            rgb = clc_st_color(angle_, ref_angle_, oblique_angle_)
            result[angle, ref_angle, :] = rgb
    return oblique_angle, result


def main_Multilayer_arrange():
    resultLUT = np.zeros((LUT_SIZE, LUT_SIZE, LUT_SIZE, 3))

    if USE_MULTIPROCESS:
        # 並列化のためのプロセスプールを作成
        with mp.Pool(processes=mp.cpu_count()) as pool:
            # oblique_angleごとに並列処理を実行
            results = pool.starmap(process_oblique_angle,
                                   [(oblique_angle,) for oblique_angle in range(LUT_SIZE)])

        # 結果を統合
        for oblique_angle, partial_result in results:
            resultLUT[oblique_angle, :, :, :] = partial_result
    else:
        for oblique_angle in range(LUT_SIZE):
            oblique_angle, result = process_oblique_angle(oblique_angle)
            resultLUT[oblique_angle, :, :, :] = result

    return resultLUT


def save_result(result):
    if SAVE_NDARRAY:
        np.save(f"../LUT/Multilayer_para3_{LUT_SIZE}", result)
    if SAVE_CSV:
        save_result_to_csv(result, f"../LUT/Multilayer_para3_{LUT_SIZE}_{PROJECT_NAME}.csv")


if __name__ == '__main__':
    print_(datetime.datetime.now(), PRINTLOG_PATH)
    start = time.perf_counter()
    resultLUT = main_Multilayer_arrange()
    end = time.perf_counter()
    save_result(resultLUT)
    print_("実行時間：" + str(end - start), PRINTLOG_PATH)



