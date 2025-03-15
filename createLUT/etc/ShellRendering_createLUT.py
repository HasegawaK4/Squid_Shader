import numpy as np
from scipy import integrate
from scipy.signal import convolve2d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cv2
#import colour
import csv
#import colour
import time
from PIL import Image
lmdmin=390
class Csv_Color:
 # C S V ファイルを読み込み, 配列に格納する
  def csv_read(self, file):
    csvfile=open(file,'r',encoding='utf-8')
    reader=csv.reader(csvfile)
    cmf = []
    cmf2 = []
    #for row in reader:
    #  print(row)
    for row in reader:
      cmf.append(row)
    csvfile.close()
    for row in cmf:
      cmf2.append([float(n) if n!='\ufeff390' else lmdmin for n in row])
    return np.array(cmf2)

cc =Csv_Color()
cmf = cc.csv_read('lin2012xyz2e_5_7sf.csv')
light_col=cc.csv_read('D65.csv')
wavelengths = cmf[:, 0]

#斜め薄膜干渉を計算しているクラス
class Thinfilm_Diag:#入射角と反射角と膜厚と屈折率と色空間と色温度
    def __init__(self, angle,d,refangle,diag_angle):#コンストラクタ
        self.n1 = 1.0 # 1層目の屈折率
        self.n3 = 1.0 # 2層目の屈折率
        self.n2=1.53
        self.diag_angle=diag_angle
        self.angle = float(angle)+self.diag_angle#傾斜分足す
        self.d=float(d)
        self.refangle=np.abs(float(refangle)-self.diag_angle)#傾斜分引く
        self.light_col = light_col
        #self.rgbspace=rgbspace
        self.clc_st_color() #構造色を計算する

    def clc_st_color(self):
        #print(self.cmf)
        wavelengths = cmf[:, 0]
        xyz = np.zeros(3)
        I = integrate.simps(cmf[:, 2], wavelengths)
        # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
        xyz[0] = integrate.simps(self.f(wavelengths, 1), wavelengths)
        xyz[1] = integrate.simps(self.f(wavelengths, 2), wavelengths)
        xyz[2] = integrate.simps(self.f(wavelengths, 3), wavelengths)
        kk = integrate.simps(self.f2(wavelengths, 2), wavelengths)
        k = 1 / kk
        xyz *= k
        self.xyz_to_widegamut(xyz)
    def xyz_to_widegamut(self, xyz):
  #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix=np.array([[1.4625,-0.1845,-0.2734],[-0.5228,1.4479,0.0681],[0.0346,-0.0958,1.2875]])
        rgb=np.dot(trans_matrix,xyz)
        
        for i in range(3):
          if rgb[i]<=0.0031308:
            rgb[i]*=12.92
          elif rgb[i]>0.0031308:
            rgb[i]=1.055*np.power(rgb[i],1/2.4)-0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def f(self, wavelengths, index):
        #r = self.BRDF(wavelengths) # 薄膜干渉の反射率を求める
        r=self.clc_reflectance(wavelengths)
        #print(self.cmf)
        col = r * cmf[:,index] * self.light_col[:,1]
        return col

    def f2(self, wavelengths, index):
        col = cmf[:,index] * self.light_col[:,1]
        return col
    
    def clc_reflectance(self,wavelength):
        rad = np.radians(self.angle)
        refrad=np.radians(self.refangle)
        cos1 = np.cos(rad)#入射角コサイン
        sin01=np.sin(rad)#入射角サイン
        sin2 = np.sin(rad)*self.n1 / self.n2
        rad2 = np.arcsin(sin2)#入射角の屈折角
        cos2 = np.cos(rad2)
        cos3 = cos1
        tan1=np.tan(rad2)#入射角の屈折角のタンジェント
        cos0=np.cos(refrad)#反射角コサイン
        sin00=np.sin(refrad)#反射角サイン
        sin0 = np.sin(refrad) / self.n2
        rad0 = np.arcsin(sin0)#反射角の屈折角
        cos4=np.cos(rad0)#反射角の屈折角のコサイン
        tan2=np.tan(rad0)#反射角の屈折角のタンジェント
# フレネル係数
        rs12 = (self.n1 * cos1 - self.n2 * cos2) / (self.n1 * cos1 +self.n2 * cos2)
        rp12 = (self.n2 * cos1 - self.n1 * cos2) / (self.n2 * cos1 + self.n1 * cos2)
        rs23 = (self.n2 * cos2 - self.n3 * cos3) / (self.n2 * cos2 + self.n3 * cos3)
        rp23 = (self.n3 * cos2 - self.n2 * cos3) / (self.n3 * cos2 + self.n2 * cos3)
        r12=(rs12+rp12)*0.5
        r23=(rs23+rp23)*0.5
        r21=1-r12#todo:これであってる？
        t12=1-r12#todo:これであってる？
        t21=r12#todo:これであってる？
        t23=1-r23#todo:これであってる？

        
        delta=self.n2*self.d*(1/cos2+1/cos4)-self.n1*self.d*(tan1+tan2)*sin01#異角度光路差
        phi = (2.0 * np.pi * delta) / wavelength # 位相差
        if self.angle<90:
            rr=np.square(np.abs(r12+t12*t21*r23*np.exp(1j*phi)/(1-r23*r21*np.exp(1j*phi))))
            return rr
        else:#self.angle>=90のとき、透過率を薄膜干渉反射率にしている
            rr=t12*t23*np.exp(1j*phi/2)/(1-r23*r12*np.exp(1j*phi))
            return rr

#斜めじゃない多層膜干渉を計算しているクラス↓
class Multilayer:#入射角と反射角と膜厚と屈折率と色空間と色温度
    def __init__(self, angle,refangle,filmnum):#コンストラクタ
        self.n_air=1
        #self.n_thick = 1.53 # アラゴナイト(炭酸カルシウム)1.53~1.69
        #self.n_thin = 1.5 # タンパク質、屈折率不明だったが高屈折率というのと、「組織透明化試薬を用いた3D 蛍光イメージングのすすめ」に1.5と記述があるため1.5
        self.n_thick=1.5
        self.n_thin=2.42

        self.angle = float(angle)
        #self.d_thick=500#アラゴナイトクリスタル層(単位:nm)
        #self.d_thin=25#タンパク質(コンキオリン)層(単位:nm)
        self.d_thick=400#アラゴナイトクリスタル層(単位:nm)
        self.d_thin=20#タンパク質(コンキオリン)層(単位:nm)
        self.a=np.pi/4#貝殻の傾斜角
        self.d=[]
        self.refangle=float(refangle)
        self.filmnum=filmnum
        self.r=[]
        self.n=[]
        self.phi=[]

        self.clc_st_color() #構造色を計算する


    def clc_st_color(self):
        #print(self.cmf)
        wavelengths = cmf[:, 0]
        xyz = np.zeros(3)
        I = integrate.simps(cmf[:, 2], wavelengths)
        # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
        self.init_param(wavelengths)
        xyz[0] = integrate.simps(self.f(wavelengths, 1), wavelengths)
        xyz[1] = integrate.simps(self.f(wavelengths, 2), wavelengths)
        xyz[2] = integrate.simps(self.f(wavelengths, 3), wavelengths)
        kk = integrate.simps(self.f2(wavelengths, 2), wavelengths)
        k = 1 / kk
        xyz *= k
        self.xyz_to_widegamut(xyz)
        #print(xyz)
    def xyz_to_srgb(self, xyz):
  #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix=np.array([[3.240970,-1.537383,-0.498611],[-0.969244,1.875968,0.041555],[0.055630,-0.203977,1.056972]])
        rgb=np.dot(trans_matrix,xyz)

        for i in range(3):
          if rgb[i]<=0.0031308:
            rgb[i]*=12.92
          elif rgb[i]>0.0031308:
            rgb[i]=1.055*np.power(rgb[i],1/2.4)-0.055

    
        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255
    def xyz_to_widegamut(self, xyz):
  #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix=np.array([[1.4625,-0.1845,-0.2734],[-0.5228,1.4479,0.0681],[0.0346,-0.0958,1.2875]])
        rgb=np.dot(trans_matrix,xyz)
        
        for i in range(3):
          if rgb[i]<=0.0031308:
            rgb[i]*=12.92
          elif rgb[i]>0.0031308:
            rgb[i]=1.055*np.power(rgb[i],1/2.4)-0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255
        #print(self.rgb)
    
    def init_param(self,wavelengths):
        self.calc_n()
        #print(self.n)
        self.calc_d()
        self.calc_r()
        self.calc_phi(wavelengths)

    def f(self, wavelengths, index):
        #r = self.BRDF(wavelengths) # 薄膜干渉の反射率を求める
        r=self.clc_final_reflectance(wavelengths)
        
        #print(r)
        #print(self.cmf)
        col = r * cmf[:,index] * light_col[:,1]
        return col

    def f2(self, wavelengths, index):
        col = cmf[:,index] * light_col[:,1]
        return col

    def calc_n(self):
        self.n.append(self.n_air)
        for i in range(self.filmnum+1):#最後の層までいれるため+1した回数ループ
            if i%2==0:
                self.n.append(self.n_thick)
            else:
                self.n.append(self.n_thin)

    def calc_d(self):
        for i in range(self.filmnum):#最後の層までいれるため+1した回数ループ
            if i%2==0:
                self.d.append(self.d_thick)
            else:
                self.d.append(self.d_thin)
        
    def calc_r(self):
        tmpangle=np.radians(self.angle)
        #tmprefangle=np.radians(self.refangle)
        tmprefracangle=np.arcsin(np.sin(tmpangle)*self.n[0]/self.n[1])
        for i in range(self.filmnum+1):#境界の数はfilmnum+1
            r_s=(self.n[i]*np.cos(tmpangle)-self.n[i+1]*np.cos(tmprefracangle))/(self.n[i]*np.cos(tmpangle)+self.n[i+1]*np.cos(tmprefracangle))
            r_p=(self.n[i+1]*np.cos(tmpangle)-self.n[i]*np.cos(tmprefracangle))/(self.n[i+1]*np.cos(tmpangle)+self.n[i]*np.cos(tmprefracangle))
            r=(r_s+r_p)/2
            self.r.append(r)
            if i!=self.filmnum:
                #print(tmpangle)
                tmpangle=tmprefracangle
                tmprefracangle=np.arcsin(self.n[i+1]*np.sin(tmprefracangle)/self.n[i+2])
                #print(tmpangle)
                #print(tmprefracangle)

    def calc_path_diff(self,n1,n2,d,angle,refrac_angle,refrac_ref_angle):#角度は弧度法
        return n2*d*(1/np.cos(refrac_angle)+1/np.cos(refrac_ref_angle))-n1*d*(np.tan(refrac_angle)+np.tan(refrac_ref_angle))*np.sin(angle)#異角度光路差

    def calc_phi(self,wavelength):#air
        tmpangle=np.radians(self.angle)
        tmp_refrac_angle=np.arcsin(self.n[0]*np.sin(tmpangle)/self.n[1])
        tmp_refrac_ref_angle=np.arcsin(self.n[0]*np.sin(np.radians(self.refangle))/self.n[1])
        for i in range(self.filmnum):
            delta=self.calc_path_diff(self.n[i],self.n[i+1],self.d[i],tmpangle,tmp_refrac_angle,tmp_refrac_ref_angle)
            phi=2*np.pi*delta/wavelength
            self.phi.append(phi)
            tmpangle=tmp_refrac_angle
            tmp_refrac_angle=np.arcsin(self.n[i+1]*np.sin(tmp_refrac_angle)/self.n[i+2])
            tmp_refrac_ref_angle=np.arcsin(self.n[i+1]*np.sin(tmp_refrac_ref_angle)/self.n[i+2])
            


        
        #if self.refangle==self.angle:
        #    print(self.phi)

    def clc_final_reflectance(self,wavelength):#層が複数になった最終的な反射率を求める
        #print(wavelength)
        # print("a")
        Rlist=[self.r[self.filmnum]]#gamma_1から順に格納
        
        #if self.angle==52.857142857142854 and self.refangle==34.285714285714285:
        #print(self.r)
        #print(self.angle,self.refangle)
        #print(self.filmnum)
        for i in range(self.filmnum,0,-1):
            gamma=(self.r[i-1]+Rlist[self.filmnum-i]*(np.cos(2*self.phi[i-1])+1j*np.sin(2*self.phi[i-1])))/(1+self.r[i-1]*Rlist[self.filmnum-i]*(np.cos(2*self.phi[i-1])+1j*np.sin(2*self.phi[i-1])))
            #print(gamma)
            Rlist.append(gamma)
            #if self.refangle==self.angle and i==1:
            #    print(self.r[i-1]+Rlist[self.filmnum-i]*(np.cos(2*self.phi[i-1])+1j*np.sin(2*self.phi[i-1])))
        #print(len(Rlist))
        #print(abs(Rlist[self.filmnum])**2)
        
        return abs(Rlist[self.filmnum])**2
    
    def clc_diagonal_reflectance(self,wavelength):#斜め多層膜の場合の反射率(とりあえず角度だけ考慮)
        theta_refrac=np.arcsin(self.n_air*np.sin(self.angle)/self.n_thick)+self.a#屈折角+貝殻の傾斜角
        

        #多層膜の計算




        """
        for i in range(self.filmnum):
            total+=tlist[i*2+0]*tlist[i*2+1]*self.clc_reflectance(nlist[i],nlist[i+1],nlist[i+2])
            """

#斜め多層膜干渉を計算するプログラム↓
class Multilayer_diag:#入射角と反射角と膜厚と屈折率と色空間と色温度
    def __init__(self, angle,refangle,filmnum,diag_angle,n_thin):#コンストラクタ, 角度は度数法で入力
        self.n_air=1
        self.n_thick = 1.53 # アラゴナイト(炭酸カルシウム)1.53~1.69
        self.n_thin = n_thin # タンパク質、屈折率不明だったが高屈折率というのと、「組織透明化試薬を用いた3D 蛍光イメージングのすすめ」に1.5と記述があるため1.7
        #self.n_thick=1.5
        #self.n_thin=2.42
        self.diag_angle=diag_angle
        if float(angle)+self.diag_angle<90:
            self.angle = float(angle)+self.diag_angle
        else:
            self.angle = 180-float(angle)-self.diag_angle
        self.d_thick=500#アラゴナイトクリスタル層(単位:nm)
        self.d_thin=25#タンパク質(コンキオリン)層(単位:nm)
        #self.d_thick=400#アラゴナイトクリスタル層(単位:nm)
        #self.d_thin=20#タンパク質(コンキオリン)層(単位:nm)
        #self.a=np.pi/4#貝殻の傾斜角
        self.d=[]
        self.refangle=np.abs(float(refangle)-self.diag_angle)
        self.filmnum=filmnum
        self.r=[]#境界の反射率のリスト
        self.t=[]#境界の透過率のリスト
        self.n=[]
        self.phi=[]

        self.clc_st_color() #構造色を計算する


    def clc_st_color(self):
        #print(self.cmf)
        wavelengths = cmf[:, 0]
        xyz = np.zeros(3)
        I = integrate.simps(cmf[:, 2], wavelengths)
        self.init_param(wavelengths)
        # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
        #if self.angle-self.diag_angle==45 and self.diag_angle+self.diag_angle==90:
        #print(self.f_intensity(wavelengths))
        xyz[0] = integrate.simps(self.f(wavelengths, 1), wavelengths)
        xyz[1] = integrate.simps(self.f(wavelengths, 2), wavelengths)
        xyz[2] = integrate.simps(self.f(wavelengths, 3), wavelengths)
        kk = integrate.simps(self.f2(wavelengths, 2), wavelengths)
        k = 1 / kk
        xyz *= k
        self.xyz_to_widegamut(xyz)

    def xyz_to_srgb(self, xyz):
  #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix=np.array([[3.240970,-1.537383,-0.498611],[-0.969244,1.875968,0.041555],[0.055630,-0.203977,1.056972]])
        rgb=np.dot(trans_matrix,xyz)

        for i in range(3):
          if rgb[i]<=0.0031308:
            rgb[i]*=12.92
          elif rgb[i]>0.0031308:
            rgb[i]=1.055*np.power(rgb[i],1/2.4)-0.055

    
        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255
    def xyz_to_widegamut(self, xyz):
  #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix=np.array([[1.4625,-0.1845,-0.2734],[-0.5228,1.4479,0.0681],[0.0346,-0.0958,1.2875]])
        rgb=np.dot(trans_matrix,xyz)
        
        for i in range(3):
          if rgb[i]<=0.0031308:
            rgb[i]*=12.92
          elif rgb[i]>0.0031308:
            rgb[i]=1.055*np.power(rgb[i],1/2.4)-0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def init_param(self,wavelengths):
        self.calc_n()
        #print(self.n)
        self.calc_d()
        self.calc_r()
        self.calc_t()
        self.calc_phi(wavelengths)
    def f(self, wavelengths, index):
        r=self.clc_final_reflectance(wavelengths)
        #print(self.cmf)

        col = r * cmf[:,index] * light_col[:,1]
        #col=r
        return col
    def f_intensity(self, wavelengths):
        r=self.clc_final_reflectance(wavelengths)
        #print(self.cmf)

        #col = r * cmf[:,index] * light_col[:,1]
        col=r
        return col

    def f2(self, wavelengths, index):
        col = cmf[:,index] * light_col[:,1]
        return col

    def calc_n(self):
        self.n.append(self.n_air)
        for i in range(self.filmnum+1):#最後の層までいれるため+1した回数ループ
            if i%2==0:
                self.n.append(self.n_thick)
            else:
                self.n.append(self.n_thin)

    def calc_d(self):
        for i in range(self.filmnum):#最後の層までいれるため+1した回数ループ
            if i%2==0:
                self.d.append(self.d_thick)
            else:
                self.d.append(self.d_thin)
        
    def calc_r(self):
        tmpangle=np.radians(self.angle)
        #tmprefangle=np.radians(self.refangle)
        tmprefracangle=np.arcsin(np.sin(tmpangle)*self.n[0]/self.n[1])
        for i in range(self.filmnum+1):#境界の数はfilmnum+1
            r_s=(self.n[i]*np.cos(tmpangle)-self.n[i+1]*np.cos(tmprefracangle))/(self.n[i]*np.cos(tmpangle)+self.n[i+1]*np.cos(tmprefracangle))
            r_p=(self.n[i+1]*np.cos(tmpangle)-self.n[i]*np.cos(tmprefracangle))/(self.n[i+1]*np.cos(tmpangle)+self.n[i]*np.cos(tmprefracangle))
            r=(r_s+r_p)/2
            self.r.append(r)
            if i!=self.filmnum:
                tmpangle=tmprefracangle
                tmprefracangle=np.arcsin(self.n[i+1]*np.sin(tmprefracangle)/self.n[i+2])
    def calc_t(self):
        tmpangle=np.radians(self.angle)
        #tmprefangle=np.radians(self.refangle)
        tmprefracangle=np.arcsin(np.sin(tmpangle)*self.n[0]/self.n[1])
        for i in range(self.filmnum+1):#境界の数はfilmnum+1
            t_s=2*self.n[i]*np.cos(tmpangle)/(self.n[i]*np.cos(tmpangle)+self.n[i+1]*np.cos(tmprefracangle))
            t_p=2*self.n[i]*np.cos(tmpangle)/(self.n[i]*np.cos(tmprefracangle)+self.n[i+1]*np.cos(tmpangle))
            t=(t_s+t_p)/2
            self.t.append(t)
            if i!=self.filmnum:
                tmpangle=tmprefracangle
                tmprefracangle=np.arcsin(self.n[i+1]*np.sin(tmprefracangle)/self.n[i+2])

    def calc_path_diff(self,n1,n2,d,angle,refrac_angle,refrac_ref_angle):#角度は弧度法
        return n2*d*(1/np.cos(refrac_angle)+1/np.cos(refrac_ref_angle))-n1*d*(np.tan(refrac_angle)+np.tan(refrac_ref_angle))*np.sin(angle)#異角度光路差

    def calc_phi(self,wavelength):#air
        tmpangle=np.radians(self.angle)
        tmp_refrac_angle=np.arcsin(self.n[0]*np.sin(tmpangle)/self.n[1])
        tmp_refrac_ref_angle=np.arcsin(self.n[0]*np.sin(np.radians(self.refangle))/self.n[1])
        for i in range(self.filmnum):
            delta=self.calc_path_diff(self.n[i],self.n[i+1],self.d[i],tmpangle,tmp_refrac_angle,tmp_refrac_ref_angle)
            phi=2*np.pi*delta/wavelength
            self.phi.append(phi)
            tmpangle=tmp_refrac_angle
            tmp_refrac_angle=np.arcsin(self.n[i+1]*np.sin(tmp_refrac_angle)/self.n[i+2])
            tmp_refrac_ref_angle=np.arcsin(self.n[i+1]*np.sin(tmp_refrac_ref_angle)/self.n[i+2])


        
        #if self.refangle==self.angle:
        #    print(self.phi)

    def clc_final_reflectance(self,wavelength):#層が複数になった最終的な反射率を求める
        if self.angle<90:
            Rlist=[self.r[self.filmnum]]#gamma_1から順に格納
            #print(self.filmnum)
            for i in range(self.filmnum,0,-1):
                gamma=(self.r[i-1]+Rlist[self.filmnum-i]*(np.cos(2*self.phi[i-1])+1j*np.sin(2*self.phi[i-1])))/(1+self.r[i-1]*Rlist[self.filmnum-i]*(np.cos(2*self.phi[i-1])+1j*np.sin(2*self.phi[i-1])))
                Rlist.append(gamma)
            return abs(Rlist[self.filmnum])**2
        else:
            #print(self.t)
            Tlist=[self.t[self.filmnum]]
            for i in range(self.filmnum,0,-1):
                gamma=(self.t[i-1]*Tlist[self.filmnum-i]*(np.cos(self.phi[i-1])+1j*np.sin(self.phi[i-1])))/(1+self.r[i-1]*Tlist[self.filmnum-i]*(np.cos(2*self.phi[i-1])+1j*np.sin(2*self.phi[i-1])))
                Tlist.append(gamma)
            return abs(Tlist[self.filmnum])**2
    
    def clc_intensity(self,filmnum):
        res=np.zeros(780-390+1)
        for lamb in range(390,781):
            sp=self.clc_final_reflectance(float(lamb))
            res[lamb-390]=sp[filmnum-1]
        return res
        

#回折を計算するプログラム↓
class Diffraction:#入射角と反射角と膜厚と屈折率と色空間と色温度
    def __init__(self, refangle,pipodnum):#コンストラクタ
        self.refangle=float(refangle)
        self.m=pipodnum#ピット数
        self.d=1#ピット幅, 単位はnm。
        self.D=8474#ピット間距離、単位はnm。
        self.DiffractionIntensity=0.02#Diffraction強度の値を調整する係数のようなもの

        self.clc_st_color() #構造色を計算する


    def clc_st_color(self):
        #print(self.cmf)
        wavelengths = cmf[:, 0]
        xyz = np.zeros(3)
        I = integrate.simps(cmf[:, 2], wavelengths)
        # wavelengthsの範囲で積分し, 反射率からx y z 表色系に変換する
        xyz[0] = integrate.simps(self.f(wavelengths, 1), wavelengths)
        xyz[1] = integrate.simps(self.f(wavelengths, 2), wavelengths)
        xyz[2] = integrate.simps(self.f(wavelengths, 3), wavelengths)
        kk = integrate.simps(self.f2(wavelengths, 2), wavelengths)
        k = 1 / kk
        xyz *= k
        self.xyz_to_widegamut(xyz)

    def xyz_to_srgb(self, xyz):
  #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix=np.array([[3.240970,-1.537383,-0.498611],[-0.969244,1.875968,0.041555],[0.055630,-0.203977,1.056972]])
        rgb=np.dot(trans_matrix,xyz)

        for i in range(3):
          if rgb[i]<=0.0031308:
            rgb[i]*=12.92
          elif rgb[i]>0.0031308:
            rgb[i]=1.055*np.power(rgb[i],1/2.4)-0.055

    
        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255
    def xyz_to_widegamut(self, xyz):
  #xyzからrgbへ変換する
        #rgb = colour.XYZ_to_sRGB(xyz) # 以下のコメントアウト内の処理を行うのと同様
        #変換行列
        trans_matrix=np.array([[1.4625,-0.1845,-0.2734],[-0.5228,1.4479,0.0681],[0.0346,-0.0958,1.2875]])
        rgb=np.dot(trans_matrix,xyz)
        
        for i in range(3):
          if rgb[i]<=0.0031308:
            rgb[i]*=12.92
          elif rgb[i]>0.0031308:
            rgb[i]=1.055*np.power(rgb[i],1/2.4)-0.055

        rgb[rgb > 1] = 1
        rgb[rgb < 0] = 0
        self.rgb = rgb * 255

    def f(self, wavelengths, index):
        #r = self.BRDF(wavelengths) # 薄膜干渉の反射率を求める
        r=self.clc_final_reflectance(wavelengths)
        #print(self.cmf)
        col = r * cmf[:,index] * light_col[:,1]
        return col
    def f_intensity(self, wavelengths):
        #r = self.BRDF(wavelengths) # 薄膜干渉の反射率を求める
        r=self.clc_final_reflectance(wavelengths)
        #print(self.cmf)
        col = r
        return col

    def f2(self, wavelengths, index):
        col = cmf[:,index] * light_col[:,1]
        return col
    
    def clc_final_reflectance(self,wavelength):
        delta_d=np.pi*self.d*np.sin(np.radians(self.refangle))/wavelength
        delta_D=np.pi*self.D*np.sin(np.radians(self.refangle))/wavelength

        I=(np.sin(delta_d)**2/delta_d**2)*(np.sin(self.m*delta_D)**2/np.sin(delta_D)**2)
        I*=self.DiffractionIntensity
        
        return I

def main_mult():
    size=64
    angle_max=90
    filmnum=1000#膜数
    imageArray=np.zeros((size,size,3),np.uint8)
    for angle in range(size):
        print(angle)
        for refangle in range(size):
            angle_Freqmethod=angle/(size-1)*angle_max
            refangle_Freqmethod=refangle/(size-1)*angle_max
            color=Multilayer(angle,refangle,filmnum)
            imageArray[angle, refangle] = [color.rgb[2], color.rgb[1], color.rgb[0]]
    cv2.imwrite('Multilayer_1000layer.jpg',imageArray)
"""
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

#回折のLUTを出力するメイン関数. 反射角, 
def main_diffraction():
    size=1024#LUTのサイズ
    angle_max=90#角度の最大値, 度数法
    imageArray=np.zeros((size,size,3),np.uint8)#空の配列
    for refangle in range(size):#反射角
        for pitnum in range(size):#ピット数
            refangle_Freqmethod=refangle/(size-1)*angle_max#0~sizeの間に正規化
            color=Diffraction(refangle,pitnum+1)#回折による色を計算
            imageArray[refangle,pitnum]=[color.rgb[2], color.rgb[1], color.rgb[0]]#Opencvで出力しているため, 順番はBGR
    cv2.imwrite('DiffractionLUT_size1024.jpg',imageArray)


#横軸を波長, 縦軸を強度とした, 回折のスペクトルをプロットするプログラム. 反射角の値が0度から90度までの91枚のグラフが出力されるので要注意
def main_Diffrac_spectrum():
    size=91#度数法
    angle_max=90#度数法
    x=np.linspace(390,780,780-390+1)
    filmnum=32#層数
    for refangle in range(size):
        refangle_=refangle/(size-1)*angle_max

        #強度を計算
        color_diffrac=Diffraction(refangle,filmnum)
        ar_diffrac=color_diffrac.f_intensity(wavelengths)
        
        #プロット
        fig, ax=plt.subplots()
        ax.plot(x,ar_diffrac,label="diffraction")
        ax.set_xlabel('wavelength(nm)')
        ax.set_ylabel('intensity')
        ax.set_title('Diffraction intensity')

        ax.legend()
        plt.savefig('spectol_diffrac_filmnum32_'+str(refangle)+'degree.png')
        plt.clf()


#斜めじゃない多層膜干渉のルックアップテーブルを出力
def main_Multilayer_3D():
    size=64
    angle_max=90
    st=''
    for refangle in range(size):#反射角
        for angle in range(size):#入射角
            for filmnum in range(size):#層数
                angle_=angle/(size-1)*angle_max#0~sizeの間に正規化
                refangle_=refangle/(size-1)*angle_max#0~sizeの間に正規化
                color=Multilayer(angle,refangle,filmnum)#多層膜干渉による色を計算
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]

                #csvファイルに書き込む文字列ファイルstに, 書き足していく
                for c in range(3):
                    st+=str(Col[c])
                    if filmnum==size and c==2:
                        pass
                    else:
                        st+=","
            st+="\n"
    with open("C:\\Users\\k2525\\Research\\MasterThesis\\Programs\\multilayer\\Multilayer.csv",mode='w') as file:#csvにルックアップテーブルを出力
        file.write(st)

#斜めの多層膜干渉のルックアップテーブルを出力
def main_Multilayer_diag():
    size=64
    angle_max=90#度数法
    diag_angle=3.335
    n_thin=1.25#コンキオリン層の屈折率. ただ, コンキオリン層の屈折率について書いてある書籍がないため, 
    st=''
    for refangle in range(size):
        for angle in range(size):
            print(refangle,angle)
            for filmnum in range(size):
                angle_=angle/(size-1)*angle_max
                refangle_=refangle/(size-1)*angle_max
                #print(refangle_,angle_)
                color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle,n_thin)
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]

                #csvファイルに書き込む文字列ファイルstに, 書き足していく
                for c in range(3):
                    st+=str(Col[c])
                    if filmnum==size and c==2:
                        pass
                    else:
                        st+=","
            st+="\n"
    with open("C:\\Users\\k2525\\Research\\MasterThesis\\Programs\\multilayer\\Multilayer.csv",mode='w') as file:
        file.write(st)

#入射角(=反射角), 層数, 傾斜角をパラメータとした, 多層膜干渉のルックアップテーブルを出力
def main_Multilayer_diag_para2():#入射角(=反射角), 層数, 傾斜角をパラメータとする
    size=64
    angle_max=90#度数法
    #diag_angle=0#度数法
    #diag_angle=3.335
    n_thin=1.25
    st=''
    for diag_angle in range(size):
        for angle in range(size):
            for filmnum in range(size):
                angle_=angle/(size-1)*angle_max
                diag_angle_=diag_angle/(size-1)*angle_max
                refangle_=angle_
                color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle_,n_thin)
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
                for c in range(3):
                    st+=str(Col[c])
                    if filmnum==size and c==2:
                        pass
                    else:
                        st+=","
            st+="\n"
    with open("C:\\Users\\k2525\\Research\\MasterThesis\\Programs\\multilayer\\Multilayer_para2.csv",mode='w') as file:
        file.write(st)
#入射角, 反射角, 傾斜角(or層数)をパラメータとした, 多層膜干渉のルックアップテーブルを出力
def main_Multilayer_diag_para3():
    size=64
    angle_max=90#度数法
    #diag_angle=0#度数法
    #diag_angle=3.335
    filmnum=32
    n_thin=1.25
    st=''
    for diag_angle in range(size):
        for angle in range(size):
            print(diag_angle,angle)
            for refangle in range(size):
                angle_=angle/(size-1)*angle_max
                diag_angle_=diag_angle/(size-1)*angle_max
                refangle_=refangle/(size-1)*angle_max
                #print(refangle_,angle_)
                color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle_,n_thin)
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
                for c in range(3):
                    st+=str(Col[c])
                    if filmnum==size and c==2:
                        pass
                    else:
                        st+=","
            st+="\n"
    with open("C:\\Users\\k2525\\Research\\MasterThesis\\Programs\\multilayer\\Multilayer_para3.csv",mode='w') as file:
        file.write(st)

#斜めの多層膜干渉の2時限ルックアップテーブルを出力している. 反射角と入射角をパラメータとしている
def main_Multilayer_diag_2D():
    size=64
    angle_max=90#度数法
    diag_angle=0#度数法
    st=''
    filmnum=1
    imageArray=np.zeros((size,size,3),np.uint8)
    for refangle in range(size):
        for angle in range(size):
            angle_=angle/(size-1)*angle_max
            refangle_=refangle/(size-1)*angle_max
            color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle)
            Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
            imageArray[refangle,angle]=[color.rgb[2], color.rgb[1], color.rgb[0]]
    cv2.imwrite('Multilayer_diag.jpg',imageArray)

#斜めの多層膜干渉の2時限ルックアップテーブルを出力している. 反射角と入射角をパラメータとし, コンキオリン層の屈折率を1~2でに100段階で変化させた100枚が出力. 
def main_Multilayer_diag_2D_refracindex():
    size=64
    angle_max=90#度数法
    diag_angle=0#度数法
    st=''
    filmnum=1
    imageArray=np.zeros((size,size,3),np.uint8)
    for refracind in range(100):
        refrac=1+1/(refracind+1)
        print(refrac)
        for refangle in range(size):
            #print(refangle)
            for angle in range(size):
                angle_=angle/(size-1)*angle_max
                refangle_=refangle/(size-1)*angle_max
                #print(angle_,refangle_)
                color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle,refracind)
                #color=Multilayer(angle_,refangle_,filmnum)
                Col = [color.rgb[0], color.rgb[1], color.rgb[2]]
                imageArray[refangle,angle]=[color.rgb[2], color.rgb[1], color.rgb[0]]
        cv2.imwrite('2D_LUT\\Multilayer_diag_n='+str(refrac)+'.jpg',imageArray)

##回折と多層膜のスペクトルを重ねてグラフとして表示
def main_Multilayer_diag_spectrum():
    size=91#0から90度で求めたい
    angle_max=90#度数法
    diag_angle=3.335#度数法
    n_thin=1.25
    st=''
    x=np.linspace(390,780,780-390+1)
    filmnum=32
    for refangle in range(size):
        for angle in range(size):
            print(refangle,angle)
            angle_=angle/(size-1)*angle_max
            refangle_=refangle/(size-1)*angle_max
            color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle,n_thin)
            color_nodiag=Multilayer_diag(angle_,refangle_,filmnum,0,n_thin)
            color_diffrac=Diffraction(refangle,filmnum)
            ar=color.f_intensity(wavelengths)
            ar_nodiag=color_nodiag.f_intensity(wavelengths)
            ar_diffrac=color_diffrac.f_intensity(wavelengths)
            
            fig, ax=plt.subplots()
            ax.plot(x,ar,label='diagonal multifilm')
            ax.plot(x,ar_nodiag,label='multifilm')
            ax.plot(x,ar_diffrac,label="diffraction")
            ax.set_xlabel('wavelength(nm)')
            ax.set_ylabel('intensity')
            ax.set_title('stcol intensity')

            ax.legend()
            #plt.show()
            plt.savefig('spectol_filmnum32_'+str(angle)+'degree_'+str(refangle)+'degree.png')
            plt.clf()
    

    # refangle=90
    # angle=45
    # angle_=angle/(size-1)*angle_max
    # refangle_=refangle/(size-1)*angle_max
    # #color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle,n_thin)
    # color_nodiag=Multilayer_diag(angle_,refangle_,filmnum,0,n_thin)
    # #ar=color.clc_intensity(filmnum)
    # ar_nodiag=color_nodiag.f_intensity(wavelengths)
    # #print(ar_nodiag)

#シャープさんに出すカラーマップの作成
def main_Multilayer_diag_colormap():
    size=91#0から90度で求めたい

    #colormap_num=41
    colormap_num=24
    #層数の時は以下をコメントアウト外す
    #colormap_num=64
    filmnum=32
    diag_angle=3.33
    #n_thin=1.25
    nthin_min=1
    nthin_max=2

    angle_max=90#度数法
    #diag_angle=3.335#度数法
    colormap=np.zeros((91,91,3))
    for i in range(colormap_num):
        print(i)
        for refangle in range(size):
            for angle in range(size):
                #傾斜角と層数とコンキオリン屈折率で以下のコメントアウト変える
                #diag_angle=i/(colormap_num-1)*angle_max
                #filmnum=i
                n_thin=nthin_min+(nthin_max-nthin_min)/(colormap_num-1)*i

                angle_=angle/(size-1)*angle_max
                refangle_=refangle/(size-1)*angle_max
                color=Multilayer_diag(angle_,refangle_,filmnum,diag_angle,n_thin)
                colormap[refangle,angle]=np.array([color.rgb[2], color.rgb[1], color.rgb[0]])
        cv2.imwrite('colormap_konnkiorinn_width_inangle_height_refangle\\colormap_'+str(n_thin)+'.jpg',colormap) 



#以下, 実行したいメイン関数のコメントアウトを外して動かす

#main_Multilayer_3D()
#main_Multilayer_diag()
#main_Multilayer_diag_para2()
#main_Multilayer_diag_spectrum()
#main_Diffrac_spectrum()
main_Multilayer_diag_2D()
#main_mult()
#main_diffraction()
#main_Multilayer_diag_2D_refracindex()
#main_Multilayer_diag_colormap()
#main_Multilayer_diag_para3()