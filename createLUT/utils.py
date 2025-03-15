import glob
import os
import re
import csv
import numpy as np

def print_(arc_str, path):
    print(arc_str)
    with open(path, 'a') as f:
        print(arc_str, file=f)


def make_filename_by_seq(dirname, filename, seq_digit=3):
    if not os.path.exists(dirname):
        # ディレクトリが存在しない場合、ディレクトリを作成する
        os.makedirs(dirname)

    filename_without_ext, ext = os.path.splitext(filename)

    pattern = f"{filename_without_ext}_([0-9]*){ext}"
    prog = re.compile(pattern)

    files = glob.glob(
        os.path.join(dirname, f"{filename_without_ext}_[0-9]*{ext}")
    )

    max_seq = -1
    for f in files:
        m = prog.match(os.path.basename(f))
        if m:
            max_seq = max(max_seq, int(m.group(1)))

    new_filename = f"{dirname}/{filename_without_ext}_{max_seq + 1:0{seq_digit}}{ext}"

    return new_filename


def judge_value(arg, dtype, error_massage=""):
    if dtype == "int":
        return int(arg.split("=")[1])
    if dtype == "float":
        return float(arg.split("=")[1])
    if dtype == "bool":
        return judge_torf(arg, error_massage)
    if dtype == "sty":
        return arg.split("=")[1]


def judge_torf(arg, error_massage=""):
    value = arg.split("=")[1].lower()
    if value in ["true", "1"]:
        return True
    elif value in ["false", "0"]:
        return False
    else:
        raise ValueError(f"{error_massage} must be a boolean (True/False or 1/0)")


def xyz_to_widegamut(xyz):
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
    rgb = rgb * 255
    return rgb


# C S V ファイルを読み込み, 配列に格納する
def csv_read(file, lmdmin):
    csvfile = open(file, 'r', encoding='utf-8')
    reader = csv.reader(csvfile)
    cmf = []
    cmf2 = []
    #for row in reader:
    #  print(row)
    for row in reader:
        cmf.append(row)
    csvfile.close()
    for row in cmf:
        cmf2.append([float(n) if n != '\ufeff390' else lmdmin for n in row])
    return np.array(cmf2)


def save_result_to_csv(result, filename):
    size = result.shape[0]# resultLUTのサイズを取得
    st = ''

    # resultLUTの内容をst形式に変換
    for diag_angle in range(size):
        for angle in range(size):
            for refangle in range(size):
                # 各colorの値をst形式に変換して追加
                for c in range(3):
                    st += str(result[diag_angle, angle, refangle, c].item())
                    st += ","
            st += "\n"  # refangleのループが終わったら改行を追加

    # CSVファイルに書き込む
    with open(filename, mode='w') as file:
        file.write(st)