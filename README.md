# high-entropy
## 程序简介

本程序采用python3编写，具体用到的第三方库有numpy, scipy, os,yaml,csv,six。如有问题，及时联系J.Pei(J.Pei@foxmail.com)。

本程序源代码托管在github上面，如需要查看最新版本程序，请移步至: [https://github.com/13skeleton/high-entropy](https://github.com/13skeleton/high-entropy)

## 程序下载方法

**Windows:**

1. 打开 [https://github.com/13skeleton/high-entropy](https://github.com/13skeleton/high-entropy) 链接，点击"clone or download"按钮，将zip文件下载至本地

2. 解压缩.zip文件，打开bin文件夹，双击运行".exe"文件。

**linux:**

> 暂时还未打包，如果想偷懒可以自己打包，如果会构建python环境，直接用python运行也行。自己解决吧，不会的，联系我。

## 程序使用方法

1. 在程序运行目录下准备一个“xxxxx.yaml”文件。

具体格式如下:

```yaml
#必填设置
Sample_Name: Bi0.3Sb1.7Te #样品名称，无实际意义，仅为区分
Component_Name: [GeTe,MnTe,SnTe] # 必填设置


##可选设置
#组元数目，如果没有此参数，默认根据Component_Name计算。
Component_Number: 3

#选填设置
Component_Proportion: [0.5,0.25,0.25] # 如果不输入该选项，默认上述组元等分。

#可选参数"Lattice_Type"，默认为cubic
Lattice_Type: cubic # cubic, orthorhombic,hexagonal，仅执行首字母判定，首字母为c/C表示cubic，首字母为o/O表示orth，首字母为h/H，表示Hexagonal；默认为cubic点阵，默认项可不填。

#可选参数，晶胞常数，单位A, 如果设定了溶解度因子，该项可不设定；如果设置了有效晶胞参数，也可不设定此项。
Lattice_Constant: [[1.08,1.06,1.09],[3,4,5],[3,5,6]] #for cubic,only need a values; for orthorhombic or hexagonal，need a、b、c value;

#可选参数，有效晶胞常数，单位A，如果设定了溶解度因子，该项可不设定；
#Effective_Lattice_Constant: [a,b,c,d,e....] #有效晶胞常数，每种物质有一个常数a


#可选参数，剪切模量，单位: GPa，如果设定了溶解度因子，该项可不设定；如果设置了平均剪切模量，也可不设定此项。
Shear_Modulus: [40,60,80]

# 平均剪切模量，单位: GPa, 可选参数，如果设定了溶解度因子，该项可不设定；
#Average_Shear_Modulus: 142.0

# 结构基元数目，可选参量，默认值为1
Z_Value: 4    # 结构基元数


#可选参数，无量纲的M值，默认值为7.34。
M_Value: 7.34


# 溶解度因子delta,单位: GPa A^3。如果设置了该参数，上述晶胞常数，剪切模量等参数不起效。
#Solubility_Factor: 1230   可直接给出溶解度因子强行计算。

#可选参数，计算某个温度的值，默认值为300K
Temperature: 400.0

#可选参数组，计算温度区间
# Temperature:
   # Start_Temperature: 300  #起始温度 单位: K
   # End_Temperature: 800    #终止温度 单位: K
   # Interval_Temperature: 1 # 温度间隔 单位: K

```

> 该输入文件遵循yaml的书写规范。可自行调整。

2. 运行程序

```bash
python high-entropy.py
```

