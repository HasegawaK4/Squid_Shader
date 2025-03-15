import numpy as np
from scipy import integrate
from scipy.signal import convolve2d
import matplotlib
import matplotlib.pyplot as plt
import cv2
#import colour
import csv
#import colour
import time
from PIL import Image
from utils import *

matplotlib.use('Agg')

lmdmin = 390

save_name = f"original"

printlog_path = make_filename_by_seq("./printlog", f"{save_name}.txt")

cmf = csv_read('lin2012xyz2e_5_7sf.csv', lmdmin)
light_col = csv_read('D65.csv', lmdmin)
wavelengths = cmf[:, 0]


#多層膜干渉を計算するクラス　引数を入射角、反射角、膜厚レベルにしている
class Multilayer:  #入射角と反射角と膜厚と屈折率と色空間と色温度
    def __init__(self, angle, refangle, para_d):  #コンストラクタ

        self.angle = float(angle)
        if self.angle == 90:  # 90度の場合
           self.angle = 89.999
        self.refangle = float(refangle)
        self.para_d = float(para_d)

        #屈折率　thinが細胞外液層、thickがタンパク質層である
        self.n_thin = 1.34
        self.n_thick = 1.44


        if abs(self.para_d) < 1e-6:
            self.para_d = 0
        elif abs(self.para_d - 1) < 1e-6:
            self.para_d = 1
            
        #膜厚(単位:nm)　thinが細胞外液層
        self.d_thin_Max = 120
        self.d_thin_min = 50
        self.d_thin_ave = 95 #生細胞95, 化学固定細胞87

        self.gamma_thin = np.log2((self.d_thin_Max - self.d_thin_min)/(self.d_thin_ave - self.d_thin_min))
        #self.gamma_thin = 1 #ガンマ補正アリのときは無効にする
        self.d_thin = (self.para_d**self.gamma_thin) * (self.d_thin_Max - self.d_thin_min) + self.d_thin_min  
        
        #thickがタンパク質層
        self.d_thick_Max = 160
        self.d_thick_min = 80
        self.d_thick_ave = 120 #生細胞120, 化学固定細胞115
        
        self.gamma_thick = np.log2((self.d_thick_Max - self.d_thick_min)/(self.d_thick_ave - self.d_thick_min))
        #self.gamma_thick = 1 #ガンマ補正アリのときは無効にする
        self.d_thick = (self.para_d**self.gamma_thick) * (self.d_thick_Max - self.d_thick_min) + self.d_thick_min  

        if self.para_d == 0:
            self.d_thin = self.d_thin_min
            self.d_thick = self.d_thick_min
        elif self.para_d == 1:
            self.d_thin = self.d_thin_Max
            self.d_thick = self.d_thick_Max

        self.d = []        
        self.filmnum = 15 #N=6→11, N=7→13, N=8→15, N=9→17
        self.r = []
        self.n = []
        self.phi = []


        #m3d5は最初のシリーズ　dL=120-50, dH=150-80
        #m3d6は生細胞、化学固定細胞の推定最大最小値にのっとった　dL=120-50, dH=160-80
        #m3d7はTEMによる推定最大最小値にのっとった　dL=100-18, dH=118-53
        #m3d8は生細胞にのっとった　dL=120-50, dH=160-80かつ平均値dL=95,dH=120を通るようガンマ補正
        #m3d9は化学固定細胞にのっとった　dL=120-50, dH=160-80かつ平均値dL=87,dH=115を通るようガンマ補正

        self.clc_st_color()  #構造色を計算する

    def clc_st_color(self):
        #print(self.cmf)
        wavelengths = cmf[:, 0]
        xyz = np.zeros(3)
        I = integrate.simpson(cmf[:, 2], x=wavelengths)
        # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
        self.init_param(wavelengths)
        xyz[0] = integrate.simpson(self.f(wavelengths, 1), x=wavelengths)
        xyz[1] = integrate.simpson(self.f(wavelengths, 2), x=wavelengths)
        xyz[2] = integrate.simpson(self.f(wavelengths, 3), x=wavelengths)
        kk = integrate.simpson(self.f2(wavelengths, 2), x=wavelengths)
        k = 1 / kk
        xyz *= k
        self.xyz_to_widegamut(xyz)
        #self.xyz_to_srgb(xyz)
        #print(xyz)

    def xyz_to_srgb(self, xyz):
        #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix = np.array(
            [[3.240970, -1.537383, -0.498611], [-0.969244, 1.875968, 0.041555], [0.055630, -0.203977, 1.056972]])
        rgb = np.dot(trans_matrix, xyz)

        for i in range(3):
            if rgb[i] <= 0.0031308:
                rgb[i] *= 12.92
            elif rgb[i] > 0.0031308:
                rgb[i] = 1.055 * np.power(rgb[i], 1 / 2.4) - 0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def xyz_to_widegamut(self, xyz):
        #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix = np.array([[1.4625, -0.1845, -0.2734], [-0.5228, 1.4479, 0.0681], [0.0346, -0.0958, 1.2875]])
        rgb = np.dot(trans_matrix, xyz)

        for i in range(3):
            if rgb[i] <= 0.0031308:
                rgb[i] *= 12.92
            elif rgb[i] > 0.0031308:
                rgb[i] = 1.055 * np.power(rgb[i], 1 / 2.4) - 0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255
        #print(self.rgb)

    def init_param(self, wavelengths):
        self.calc_n()
        #print(self.n)
        self.calc_d()
        self.calc_r()
        self.calc_phi(wavelengths)

    def f(self, wavelengths, index):
        #r = self.BRDF(wavelengths) # 薄膜干渉の反射率を求める
        r = self.clc_final_reflectance(wavelengths)#外に出して計算軽くできそう

        #print(r)
        #print(self.cmf)
        col = r * cmf[:, index] * light_col[:, 1]
        return col

    def f2(self, wavelengths, index):
        col = cmf[:, index] * light_col[:, 1]
        return col

    def calc_n(self):
        self.n.append(self.n_thin)
        for i in range(self.filmnum + 1):  #最後の層までいれるため+1した回数ループ
            if i % 2 == 0:
                self.n.append(self.n_thick)
            else:
                self.n.append(self.n_thin)
        #print(self.n)

    def calc_d(self):
        for i in range(self.filmnum):  #膜厚が必要な1～N層分
            if i % 2 == 0:
                self.d.append(self.d_thick)
            else:
                self.d.append(self.d_thin)
        #print(self.d)

    def calc_r(self):
        tmpangle = np.radians(self.angle)
        #tmprefangle=np.radians(self.refangle)
        tmprefracangle = np.arcsin(np.sin(tmpangle) * self.n[0] / self.n[1])#スネル
        for i in range(self.filmnum + 1):  #境界の数はfilmnum+1
            r_s = (self.n[i] * np.cos(tmpangle) - self.n[i + 1] * np.cos(tmprefracangle)) / (
                    self.n[i] * np.cos(tmpangle) + self.n[i + 1] * np.cos(tmprefracangle))
            r_p = (self.n[i + 1] * np.cos(tmpangle) - self.n[i] * np.cos(tmprefracangle)) / (
                    self.n[i + 1] * np.cos(tmpangle) + self.n[i] * np.cos(tmprefracangle))
            r = (r_s + r_p) / 2
            self.r.append(r)
            if i != self.filmnum:
                #print(i,tmpangle,tmprefracangle)
                tmpangle = tmprefracangle
                tmprefracangle = np.arcsin(np.sin(tmpangle) * self.n[i + 1] / self.n[i + 2])
                #print(i,tmpangle,tmprefracangle)
        #print(self.r)

    def calc_path_diff(self, n0, n1, d, angle, refrac_angle, refrac_ref_angle):  #角度は弧度法 デルタの中身
        return n1 * d * (1 / np.cos(refrac_angle) + 1 / np.cos(refrac_ref_angle)) - n0 * d * (
                np.tan(refrac_angle) + np.tan(refrac_ref_angle)) * np.sin(angle)  #異角度光路差

    def calc_phi(self, wavelength):  
        tmpangle = np.radians(self.angle)
        #tmp_refrac_angle = np.arcsin(self.n[0] * np.sin(tmpangle) / self.n[1]) #元の式
        tmp_refrac_angle = np.arcsin(np.clip(np.sin(tmpangle) * self.n[0] / self.n[1], -1.0, 1.0)) #GPT補正式
        tmp_refrac_ref_angle = np.arcsin(self.n[0] * np.sin(np.radians(self.refangle)) / self.n[1])
        for i in range(self.filmnum):
            delta = self.calc_path_diff(self.n[i], self.n[i + 1], self.d[i], tmpangle, tmp_refrac_angle,
                                        tmp_refrac_ref_angle)
            phi = 2 * np.pi * delta / wavelength
            self.phi.append(phi)
            
            tmpangle = tmp_refrac_angle
            tmp_refrac_angle = np.arcsin(self.n[i + 1] * np.sin(tmp_refrac_angle) / self.n[i + 2])
            tmp_refrac_ref_angle = np.arcsin(self.n[i + 1] * np.sin(tmp_refrac_ref_angle) / self.n[i + 2])

        #if self.refangle==self.angle:
        #    print(self.phi)

    def clc_final_reflectance(self, wavelength):  #層が複数になった最終的な反射率を求める
        #print(wavelength)
        # print("a")
        Rlist = [self.r[self.filmnum]]  #gamma_0から順に格納
        #print(Rlist)

        #if self.angle==52.857142857142854 and self.refangle==34.285714285714285:
        #print(self.r)
        #print(self.angle,self.refangle)
        #print(self.filmnum)
        for i in range(self.filmnum, 0, -1):#i=15→1層まで
            #tmp_phi = np.cos(2 * self.phi[i - 1]) + 1j * np.sin(2 * self.phi[i - 1])#*2要らなそう
            tmp_phi = np.cos(self.phi[i - 1]) + 1j * np.sin(self.phi[i - 1])#*2無しver
            gamma = (self.r[i - 1] + Rlist[self.filmnum - i] * tmp_phi) / (1 + self.r[i - 1] * Rlist[self.filmnum - i] * tmp_phi)
            #print(gamma)
            Rlist.append(gamma)
            #if self.refangle==self.angle and i==1:
            #    print(self.r[i-1]+Rlist[self.filmnum-i]*(np.cos(2*self.phi[i-1])+1j*np.sin(2*self.phi[i-1])))
        #print(len(Rlist))

        return abs(Rlist[self.filmnum]) ** 2

    def clc_diagonal_reflectance(self, wavelength):  #斜め多層膜の場合の反射率(とりあえず角度だけ考慮)
        theta_refrac = np.arcsin(self.n_thin * np.sin(self.angle) / self.n_thick) + self.a  #屈折角+貝殻の傾斜角

        #多層膜の計算

        """
        for i in range(self.filmnum):
            total+=tlist[i*2+0]*tlist[i*2+1]*self.clc_reflectance(nlist[i],nlist[i+1],nlist[i+2])
            """
"""
#基となる3次元LUT生成プログラム
def main_Multilayer_3D():
    size = 64
    angle_max = 90
    st = ''
    for refangle in range(size):  #反射角
        for angle in range(size):  #入射角
            print(refangle,angle)
            for filmnum in range(size):  #層数
                angle_ = angle / (size - 1) * angle_max  #0~sizeの間に正規化
                refangle_ = refangle / (size - 1) * angle_max  #0~sizeの間に正規化
                color = Multilayer(angle, refangle, filmnum)  #多層膜干渉による色を計算
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]

                #csvファイルに書き込む文字列ファイルstに, 書き足していく
                for c in range(3):
                    st += str(Col[c])
                    if filmnum == size and c == 2:
                        pass
                    else:
                        st += ","
            st += "\n"
    with open("C:\\Users\\admin\\Documents\\1研究室\\Unity\\22T0001-YAsahina-Thesis-master\\createLUT\\Multilayer_try.csv",
              mode='w') as file:  #csvにルックアップテーブルを出力
        file.write(st)
        """
        
#膜厚レベルに対応させた3次元LUT生成プログラム
def main_myMultilayer_3D():
    size = 64
    angle_max = 90
    st = ''
    for refangle in range(size):  #反射角
        refangleM = refangle / (size - 1) * angle_max  #0~90の間に正規化
        for angle in range(size):  #入射角
            angleM = angle / (size - 1) * angle_max    #0~90の間に正規化
            for para_d in range(size):#膜厚レベル    
                para_dM = para_d / (size - 1)          #0~1の間に正規化
                print(angle,refangle,para_d)
                color = Multilayer(angleM, refangleM, para_dM)  #多層膜干渉による色を計算
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]

                #csvファイルに書き込む文字列ファイルstに, 書き足していく
                for c in range(3):
                    st += str(Col[c])
                    if para_d == size-1 and c == 2:
                        pass
                    else:
                        st += ","
            st += "\n"
    with open("C:\\Users\\admin\\Documents\\1研究室\\Unity\\22T0001-YAsahina-Thesis-master\\createLUT\\Multilayer_m3d8-N8-new2.csv",
              mode='w') as file:  #csvにルックアップテーブルを出力
        file.write(st)
        #m3d5は最初のシリーズ　dL=120-50, dH=150-80
        #m3d6は生細胞、化学固定細胞の推定最大最小値にのっとった　dL=120-50, dH=160-80
        #m3d7はTEMによる推定最大最小値にのっとった　dL=100-18, dH=118-53
        #m3d8は生細胞にのっとった　dL=120-50, dH=160-80かつ平均値dL=95,dH=120を通るようガンマ補正
        #m3d9は化学固定細胞にのっとった　dL=120-50, dH=160-80かつ平均値dL=87,dH=115を通るようガンマ補正

"""
#基となる2次元LUT生成プログラム
def main_Multilayer_2D():
    size = 64
    angle_max = 90  #度数法
    st = ''
    filmnum = 15
    imageArray = np.zeros((size, size, 3), np.uint8)
    for refangle in range(size):
        for angle in range(size):
            print(refangle,angle)
            angle_ = angle / (size - 1) * angle_max
            refangle_ = refangle / (size - 1) * angle_max
            color = Multilayer(angle, refangle, filmnum)
            Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
            imageArray[refangle, angle] = [color.rgb[2], color.rgb[1], color.rgb[0]]
    cv2.imwrite('Multilayer_2d1.jpg', imageArray)
"""
#膜厚レベルを手入力にした二次元LUT生成プログラム
def main_myMultilayer_2D():
    size = 64
    angle_max =  90 #度数法
    para_d = 0.01 #膜厚レベル(０～１で選択)
    st = ''
    imageArray = np.zeros((size, size, 3), np.uint8)
    for refangle in range(size):
        refangleM = refangle / (size - 1) * angle_max
        for angle in range(size):
            angleM = angle / (size - 1) * angle_max            
            print(angle, refangle)
            color = Multilayer(angleM, refangleM, para_d)
            Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
            imageArray[refangle, angle] = [color.rgb[2], color.rgb[1], color.rgb[0]]
    cv2.imwrite('Multilayer_m2d8N8new_0.01.jpg', imageArray)
    print(para_d)








        

"""
#透過を考慮した内容?
class Multilayer_diag:  #入射角と反射角と膜厚と屈折率と色空間と色温度
    def __init__(self, angle, refangle, filmnum, diag_angle, n_thin):  #コンストラクタ, 角度は度数法で入力
        self.n_air = 1
        self.n_thick = 1.53  # アラゴナイト(炭酸カルシウム)1.53~1.69
        self.n_thin = n_thin  # タンパク質、屈折率不明だったが高屈折率というのと、「組織透明化試薬を用いた3D 蛍光イメージングのすすめ」に1.5と記述があるため1.7
        #self.n_thick=1.5
        #self.n_thin=2.42
        self.diag_angle = diag_angle
        if float(angle) + self.diag_angle < 90:
            self.angle = float(angle) + self.diag_angle
        else:
            self.angle = 180 - float(angle) - self.diag_angle
        self.d_thick = 500  #アラゴナイトクリスタル層(単位:nm)
        self.d_thin = 25  #タンパク質(コンキオリン)層(単位:nm)
        #self.d_thick=400#アラゴナイトクリスタル層(単位:nm)
        #self.d_thin=20#タンパク質(コンキオリン)層(単位:nm)
        #self.a=np.pi/4#貝殻の傾斜角
        self.d = []
        self.refangle = np.abs(float(refangle) - self.diag_angle)
        self.filmnum = filmnum
        self.r = []  #境界の反射率のリスト
        self.t = []  #境界の透過率のリスト
        self.n = []
        self.phi = []

        self.clc_st_color()  #構造色を計算する

    def clc_st_color(self):
        #print(self.cmf)
        wavelengths = cmf[:, 0]
        xyz = np.zeros(3)
        I = integrate.simpson(cmf[:, 2], x=wavelengths)
        self.init_param(wavelengths)
        # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
        #if self.angle-self.diag_angle==45 and self.diag_angle+self.diag_angle==90:
        #print(self.f_intensity(wavelengths))
        xyz[0] = integrate.simpson(self.f(wavelengths, 1), x=wavelengths)
        xyz[1] = integrate.simpson(self.f(wavelengths, 2), x=wavelengths)
        xyz[2] = integrate.simpson(self.f(wavelengths, 3), x=wavelengths)
        kk = integrate.simpson(self.f2(wavelengths, 2), x=wavelengths)
        k = 1 / kk
        xyz *= k
        self.xyz_to_widegamut(xyz)

    def xyz_to_srgb(self, xyz):
        #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix = np.array(
            [[3.240970, -1.537383, -0.498611], [-0.969244, 1.875968, 0.041555], [0.055630, -0.203977, 1.056972]])
        rgb = np.dot(trans_matrix, xyz)

        for i in range(3):
            if rgb[i] <= 0.0031308:
                rgb[i] *= 12.92
            elif rgb[i] > 0.0031308:
                rgb[i] = 1.055 * np.power(rgb[i], 1 / 2.4) - 0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def xyz_to_widegamut(self, xyz):
        #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix = np.array([[1.4625, -0.1845, -0.2734], [-0.5228, 1.4479, 0.0681], [0.0346, -0.0958, 1.2875]])
        rgb = np.dot(trans_matrix, xyz)

        for i in range(3):
            if rgb[i] <= 0.0031308:
                rgb[i] *= 12.92
            elif rgb[i] > 0.0031308:
                rgb[i] = 1.055 * np.power(rgb[i], 1 / 2.4) - 0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def init_param(self, wavelengths):
        self.calc_n()
        #print(self.n)
        self.calc_d()
        self.calc_r()
        self.calc_t()
        self.calc_phi(wavelengths)

    def f(self, wavelengths, index):
        r = self.clc_final_reflectance(wavelengths)
        #print(self.cmf)

        col = r * cmf[:, index] * light_col[:, 1]
        #col=r
        return col

    def f_intensity(self, wavelengths):
        r = self.clc_final_reflectance(wavelengths)
        #print(self.cmf)

        #col = r * cmf[:,index] * light_col[:,1]
        col = r
        return col

    def f2(self, wavelengths, index):
        col = cmf[:, index] * light_col[:, 1]
        return col

    def calc_n(self):
        self.n.append(self.n_air)
        for i in range(self.filmnum + 1):  #最後の層までいれるため+1した回数ループ
            if i % 2 == 0:
                self.n.append(self.n_thick)
            else:
                self.n.append(self.n_thin)

    def calc_d(self):
        for i in range(self.filmnum):  
            if i % 2 == 0:
                self.d.append(self.d_thick)
            else:
                self.d.append(self.d_thin)

    def calc_r(self):
        tmpangle = np.radians(self.angle)
        #tmprefangle=np.radians(self.refangle)
        tmprefracangle = np.arcsin(np.sin(tmpangle) * self.n[0] / self.n[1])
        for i in range(self.filmnum + 1):  #境界の数はfilmnum+1
            r_s = (self.n[i] * np.cos(tmpangle) - self.n[i + 1] * np.cos(tmprefracangle)) / (
                    self.n[i] * np.cos(tmpangle) + self.n[i + 1] * np.cos(tmprefracangle))
            r_p = (self.n[i + 1] * np.cos(tmpangle) - self.n[i] * np.cos(tmprefracangle)) / (
                    self.n[i + 1] * np.cos(tmpangle) + self.n[i] * np.cos(tmprefracangle))
            r = (r_s + r_p) / 2
            self.r.append(r)
            if i != self.filmnum:
                tmpangle = tmprefracangle
                tmprefracangle = np.arcsin(self.n[i + 1] * np.sin(tmprefracangle) / self.n[i + 2])

    def calc_t(self):
        tmpangle = np.radians(self.angle)
        #tmprefangle=np.radians(self.refangle)
        tmprefracangle = np.arcsin(np.sin(tmpangle) * self.n[0] / self.n[1])
        for i in range(self.filmnum + 1):  #境界の数はfilmnum+1
            t_s = 2 * self.n[i] * np.cos(tmpangle) / (
                    self.n[i] * np.cos(tmpangle) + self.n[i + 1] * np.cos(tmprefracangle))
            t_p = 2 * self.n[i] * np.cos(tmpangle) / (
                    self.n[i] * np.cos(tmprefracangle) + self.n[i + 1] * np.cos(tmpangle))
            t = (t_s + t_p) / 2
            self.t.append(t)
            if i != self.filmnum:
                tmpangle = tmprefracangle
                tmprefracangle = np.arcsin(self.n[i + 1] * np.sin(tmprefracangle) / self.n[i + 2])

    def calc_path_diff(self, n1, n2, d, angle, refrac_angle, refrac_ref_angle):  #角度は弧度法
        return n2 * d * (1 / np.cos(refrac_angle) + 1 / np.cos(refrac_ref_angle)) - n1 * d * (
                np.tan(refrac_angle) + np.tan(refrac_ref_angle)) * np.sin(angle)  #異角度光路差

    def calc_phi(self, wavelength):  #air
        tmpangle = np.radians(self.angle)
        tmp_refrac_angle = np.arcsin(self.n[0] * np.sin(tmpangle) / self.n[1])
        tmp_refrac_ref_angle = np.arcsin(self.n[0] * np.sin(np.radians(self.refangle)) / self.n[1])
        for i in range(self.filmnum):
            delta = self.calc_path_diff(self.n[i], self.n[i + 1], self.d[i], tmpangle, tmp_refrac_angle,
                                        tmp_refrac_ref_angle)
            phi = 2 * np.pi * delta / wavelength
            self.phi.append(phi)
            tmpangle = tmp_refrac_angle
            tmp_refrac_angle = np.arcsin(self.n[i + 1] * np.sin(tmp_refrac_angle) / self.n[i + 2])
            tmp_refrac_ref_angle = np.arcsin(self.n[i + 1] * np.sin(tmp_refrac_ref_angle) / self.n[i + 2])

        #if self.refangle==self.angle:
        #    print(self.phi)

    def clc_final_reflectance(self, wavelength):  #層が複数になった最終的な反射率を求める
        if self.angle < 90:
            Rlist = [self.r[self.filmnum]]  #gamma_1から順に格納
            #print(self.filmnum)
            for i in range(self.filmnum, 0, -1):
                gamma = (self.r[i - 1] + Rlist[self.filmnum - i] * (
                        np.cos(2 * self.phi[i - 1]) + 1j * np.sin(2 * self.phi[i - 1]))) / (
                                1 + self.r[i - 1] * Rlist[self.filmnum - i] * (
                                np.cos(2 * self.phi[i - 1]) + 1j * np.sin(2 * self.phi[i - 1])))
                Rlist.append(gamma)
            return abs(Rlist[self.filmnum]) ** 2
        else:
            #print(self.t)
            Tlist = [self.t[self.filmnum]]
            for i in range(self.filmnum, 0, -1):
                gamma = (self.t[i - 1] * Tlist[self.filmnum - i] * (
                        np.cos(self.phi[i - 1]) + 1j * np.sin(self.phi[i - 1]))) / (
                                1 + self.r[i - 1] * Tlist[self.filmnum - i] * (
                                np.cos(2 * self.phi[i - 1]) + 1j * np.sin(2 * self.phi[i - 1])))
                Tlist.append(gamma)
            return abs(Tlist[self.filmnum]) ** 2

    def clc_intensity(self, filmnum):
        res = np.zeros(780 - 390 + 1)
        for lamb in range(390, 781):
            sp = self.clc_final_reflectance(float(lamb))
            res[lamb - 390] = sp[filmnum - 1]
        return res


#回折を計算するプログラム↓
class Diffraction:  #入射角と反射角と膜厚と屈折率と色空間と色温度
    def __init__(self, refangle, pipodnum):  #コンストラクタ
        self.refangle = float(refangle)
        self.m = pipodnum  #ピット数
        self.d = 1  #ピット幅, 単位はnm。
        self.D = 8474  #ピット間距離、単位はnm。
        self.DiffractionIntensity = 0.02  #Diffraction強度の値を調整する係数のようなもの

        self.clc_st_color()  #構造色を計算する

    def clc_st_color(self):
        #print(self.cmf)
        wavelengths = cmf[:, 0]
        xyz = np.zeros(3)
        I = integrate.simpson(cmf[:, 2], x=wavelengths)
        # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
        xyz[0] = integrate.simpson(self.f(wavelengths, 1), x=wavelengths)
        xyz[1] = integrate.simpson(self.f(wavelengths, 2), x=wavelengths)
        xyz[2] = integrate.simpson(self.f(wavelengths, 3), x=wavelengths)
        kk = integrate.simpson(self.f2(wavelengths, 2), x=wavelengths)
        k = 1 / kk
        xyz *= k
        self.xyz_to_widegamut(xyz)

    def xyz_to_srgb(self, xyz):
        #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix = np.array(
            [[3.240970, -1.537383, -0.498611], [-0.969244, 1.875968, 0.041555], [0.055630, -0.203977, 1.056972]])
        rgb = np.dot(trans_matrix, xyz)

        for i in range(3):
            if rgb[i] <= 0.0031308:
                rgb[i] *= 12.92
            elif rgb[i] > 0.0031308:
                rgb[i] = 1.055 * np.power(rgb[i], 1 / 2.4) - 0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def xyz_to_widegamut(self, xyz):
        #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix = np.array([[1.4625, -0.1845, -0.2734], [-0.5228, 1.4479, 0.0681], [0.0346, -0.0958, 1.2875]])
        rgb = np.dot(trans_matrix, xyz)

        for i in range(3):
            if rgb[i] <= 0.0031308:
                rgb[i] *= 12.92
            elif rgb[i] > 0.0031308:
                rgb[i] = 1.055 * np.power(rgb[i], 1 / 2.4) - 0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def f(self, wavelengths, index):
        #r = self.BRDF(wavelengths) # 薄膜干渉の反射率を求める
        r = self.clc_final_reflectance(wavelengths)
        #print(self.cmf)
        col = r * cmf[:, index] * light_col[:, 1]
        return col

    def f_intensity(self, wavelengths):
        #r = self.BRDF(wavelengths) # 薄膜干渉の反射率を求める
        r = self.clc_final_reflectance(wavelengths)
        #print(self.cmf)
        col = r
        return col

    def f2(self, wavelengths, index):
        col = cmf[:, index] * light_col[:, 1]
        return col

    def clc_final_reflectance(self, wavelength):
        delta_d = np.pi * self.d * np.sin(np.radians(self.refangle)) / wavelength
        delta_D = np.pi * self.D * np.sin(np.radians(self.refangle)) / wavelength

        I = (np.sin(delta_d) ** 2 / delta_d ** 2) * (np.sin(self.m * delta_D) ** 2 / np.sin(delta_D) ** 2)
        I *= self.DiffractionIntensity

        return I
"""


"""
def main_mult():
    size = 64
    angle_max = 90
    filmnum = 1000  #膜数
    imageArray = np.zeros((size, size, 3), np.uint8)
    for angle in range(size):
        print(angle)
        for refangle in range(size):
            angle_Freqmethod = angle / (size - 1) * angle_max
            refangle_Freqmethod = refangle / (size - 1) * angle_max
            color = Multilayer(angle, refangle, filmnum)
            imageArray[angle, refangle] = [color.rgb[2], color.rgb[1], color.rgb[0]]
    cv2.imwrite('Multilayer_1000layer.jpg', imageArray)



def main_mult_1D():
    size=64
    angle_max=90
    filmnum=5#膜数
    imageArray=np.zeros((size,size,3),np.uint8)
    for angle in range(size):
        refangle=angle
        angle_Freqmethod=angle/(size-1)*angle_max
        refangle_Freqmethod=refangle/(size-1)*angle_max

        color=Multilayer(angle,refangle,filmnum)        
        imageArray[angle, :] = [color.rgb[2], color.rgb[1], color.rgb[0]]
    cv2.imwrite('Multilayer.jpg',imageArray)

"""




















"""
#入射角(=反射角), 層数, 傾斜角をパラメータとした, 多層膜干渉のルックアップテーブルを出力
def main_Multilayer_diag_para2():  #入射角(=反射角), 層数, 傾斜角をパラメータとする
    size = 64
    angle_max = 90  #度数法
    #diag_angle=0#度数法
    #diag_angle=3.335
    n_thin = 1.25
    st = ''
    for diag_angle in range(size):
        for angle in range(size):
            for filmnum in range(size):
                angle_ = angle / (size - 1) * angle_max
                diag_angle_ = diag_angle / (size - 1) * angle_max
                refangle_ = angle_
                color = Multilayer_diag(angle_, refangle_, filmnum, diag_angle_, n_thin)
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
                for c in range(3):
                    st += str(Col[c])
                    if filmnum == size and c == 2:
                        pass
                    else:
                        st += ","
            st += "\n"
    with open("C:\\Users\\k2525\\Research\\MasterThesis\\Programs\\multilayer\\Multilayer_para2.csv", mode='w') as file:
        file.write(st)


# 入射角, 反射角, 傾斜角(or層数)をパラメータとした, 多層膜干渉のルックアップテーブルを出力
def main_Multilayer_diag_para3():
    size = 16
    angle_max = 90  # 度数法
    # diag_angle=0  # 度数法
    # diag_angle=3.335
    filmnum = 32
    n_thin = 1.25
    st = ''
    resultLUT = np.zeros((size, size, size, 3))
    for diag_angle in range(size):
        for angle in range(size):
            print(diag_angle, angle)
            for refangle in range(size):
                angle_ = angle / (size - 1) * angle_max
                diag_angle_ = diag_angle / (size - 1) * angle_max
                refangle_ = refangle / (size - 1) * angle_max
                #print(refangle_,angle_)
                color = Multilayer_diag(angle_, refangle_, filmnum, diag_angle_, n_thin)
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
                # print(color.rgb.dtype)
                resultLUT[diag_angle, angle, refangle, :] = color.rgb
                for c in range(3):
                    st += str(Col[c])
                    st += ","
            st += "\n"

    with open(f"../LUT/Multilayer_para3_{size}_original.csv", mode='w') as file:
        file.write(st)
    np.save(f"../LUT/Multilayer_para3_{size}", resultLUT)





"""


    #以下, 実行したいメイン関数のコメントアウトを外して動かす

start = time.perf_counter()



#基本的に↓のどっちか動かす
main_myMultilayer_3D()
#main_myMultilayer_2D()

#main_Multilayer_3D()
#main_Multilayer_2D()

#main_Multilayer_diag()
#main_Multilayer_diag_para2()
#main_Multilayer_diag_spectrum()
#main_Diffrac_spectrum()
#main_Multilayer_diag_2D()
#main_mult()
#main_diffraction()
#main_Multilayer_diag_2D_refracindex()
#main_Multilayer_diag_colormap()
#main_Multilayer_diag_para3()

end = time.perf_counter()
print_("実行時間：" + str(end - start), printlog_path)
