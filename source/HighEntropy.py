import scipy.constants as C
import numpy as np
import os,yaml
import csv
import six

N_A = C.N_A
R = C.R

class High_Entropy(object):
    def __init__(self,name):
        self.__sample_name = name    
        self.__component_number_set = 0
        self.__component_proportion_set = 0
        self.__temperature_set = 0
        self.__lattice_type_set = 0
        self.__lattice_constant_set = 0
        self.__effective_lattice_constant_set = 0
        self.__shear_modulus_set = 0
        self.__average_shear_modulus_set = 0
        self.__Z_value_set = 0
        self.__average_effective_lattice_constant_set =0
        self.__M_value_set = 0
        self.__solubility_factor_set = 0


    @property
    def sample_name(self):
        return self.__sample_name
    
    # 读写各组元名称
    @property
    def component_name(self):
        return self.__component_name
    @component_name.setter
    def component_name(self,value):
        self.__component_name = np.array(value)

    # 计算混合组元的种类，如Sn0.5Mn0.5Te中的组元为SnTe, MnTe。
    @property
    def component_number(self):
        if self.__component_number_set == 0:
            return len(self.component_name)
        else:
            return self.__component_number
    @component_number.setter
    def component_number(self,value):
        self.__component_number = np.array(value)
        self.__component_number_set = 1

    
    # 读写各组元比例分数
    @property
    def component_proportion(self):
        if self.__component_proportion_set == 0:
            self.__component_proportion = np.ones(len(self.component_number))
        return self.__component_proportion
    @component_proportion.setter
    def component_proportion(self,value):
        self.__component_proportion = value
        self.__component_proportion_set = 1
    
    #读写物质的温度，temperature,T均表示温度,不设置，默认为300K
    @property
    def temperature(self):
        if self.__temperature_set == 0:
            self.__temperature = 300.0
        return self.__temperature
    @temperature.setter
    def temperature(self,value):
        self.__temperature = float(value)
        self.__temperature_set =1

        
    #读写基体的空间点阵类型,默认为cubic点阵
    @property
    def lattice_type(self):
        if self.__lattice_type_set == 0:
            self.__lattice_type = "cubic"
        return self.__lattice_type
    @lattice_type.setter
    def lattice_type(self,value):
        if value.capitalize()[0] =="C":
            self.__lattice_type = "cubic"
        elif value.capitalize()[0] =="O":
            self.__lattice_type = "orthorhombic"
        elif value.capitalize()[0] =="H":
            self.__lattice_type = "hexagonal"
        else:
            self.__lattice_type = "cubic"
        self.__lattice_type_set = 1

    # 读写各组元晶胞常数
    @property
    def lattice_constant(self):
        if self.__lattice_constant_set == 1:
            return self.__lattice_constant
    @lattice_constant.setter
    def lattice_constant(self,value):
        if self.lattice_type =="cubic":
            if np.array(value).ndim ==1 and (np.array(value).shape[0] == self.component_number):
                self.__lattice_constant = np.array(value).repeat(3).reshape(len(value),3)
            elif np.array(value).ndim ==2 and (np.array(value).shape[0] == self.component_number) and (np.array(value).shape[1] == 1):
                self.__lattice_constant = np.array(value).repeat(3).reshape(len(value),3)
            elif np.array(value).ndim ==2 and (np.array(value).shape[0] == self.component_number) and (np.array(value).shape[1] == 3):
                self.__lattice_constant = np.array(value)
        elif self.lattice_type =="orthorhombic":
            if np.array(value).ndim ==1 and (np.array(value).shape[0] == 3*self.component_number):
                self.__lattice_constant = np.array(value).reshape(self.component_number,3)
            elif np.array(value).ndim ==2 and (np.array(value).shape[0] == self.component_number) and (np.array(value).shape[1] == 3):
                self.__lattice_constant = np.array(value)
        elif self.lattice_type =="hexagonal":
            if np.array(value).ndim ==1 and (np.array(value).shape[0] == self.component_number*2):
                self.__lattice_constant = np.array([np.array(value).reshape(3,2).T[0],np.array(value).reshape(3,2).T[0],np.array(value).reshape(3,2).T[1]]).T
            elif np.array(value).ndim ==2 and (np.array(value).shape[0] == self.component_number) and (np.array(value).shape[1] == 2):
                self.__lattice_constant = np.array([np.array(value).T[0],np.array(value).T[0],np.array(value).T[1]]).T
            elif np.array(value).ndim ==2 and (np.array(value).shape[0] == self.component_number) and (np.array(value).shape[1] == 3):
                self.__lattice_constant = np.array(value)
        self.__lattice_constant_set = 1

        
    # 读写各组元的有效晶胞常数
    @property
    def effective_lattice_constant(self):
        if self.__effective_lattice_constant_set == 0 and self.__lattice_constant_set == 1:
            self.__effective_lattice_constant = np.sqrt(np.sum(np.power(self.lattice_constant,2)/3,axis=1))
        else:
            self.__effective_lattice_constant = np.zeros(self.component_number)
        return self.__effective_lattice_constant
    @effective_lattice_constant.setter
    def effective_lattice_constant(self,value):
        self.__effective_lattice_constant = np.array(value)
        self.__effective_lattice_constant_set =1
    
    @property
    def difference_effective_lattice_constant(self):
        if self.__effective_lattice_constant.size != 0 and np.sum(self.effective_lattice_constant)!= 0:
            self.__difference_effective_lattice_constant= np.array(self.effective_lattice_constant)-np.array(self.effective_lattice_constant)[0]
        else:
            self.__difference_effective_lattice_constant = np.zeros(self.component_number)
        return self.__difference_effective_lattice_constant


    # 读写各组元剪切模量
    @property
    def shear_modulus(self):
        if self.__shear_modulus_set == 1:
            return self.__shear_modulus
    @shear_modulus.setter
    def shear_modulus(self,value):
        self.__shear_modulus = value
        self.__shear_modulus_set =1
    
    #读写结构基元数目
    @property
    def Z_value(self):
        if self.__Z_value_set == 0:
            self.__Z_value = 1
        return self.__Z_value
    @Z_value.setter
    def Z_value(self,value):
        self.__Z_value = value
        self.__Z_value_set = 1
    
    
    #读写平均剪切模量
    @property
    def average_shear_modulus(self):
        if self.__average_shear_modulus_set == 0 and self.__shear_modulus_set == 1:
            self.__average_shear_modulus = np.average(self.shear_modulus)
        return self.__average_shear_modulus
    @average_shear_modulus.setter
    def average_shear_modulus(self,value):
        self.__average_shear_modulus = value
        self.__average_shear_modulus_set = 1
    
    #读写平均有效晶胞常数
    @property
    def average_effective_lattice_constant(self):
        if self.__average_effective_lattice_constant_set == 0:
            self.__average_effective_lattice_constant=np.average(self.effective_lattice_constant)   
        return self.__average_effective_lattice_constant

    @average_effective_lattice_constant.setter
    def average_effective_lattice_constant(self,value):
        self.__average_effective_lattice_constant = value
        self.__average_effective_lattice_constant_set = 1
    

    #读写无量纲常数M，默认值为7.34
    @property
    def M_value(self):
        if self.__M_value_set == 0:
            self.__M_value = 7.34
        return self.__M_value
    @M_value.setter
    def M_value(self,value):
        self.__M_value = value
        self.__M_value_set = 1
    
    #计算溶解度因子delta(单位:GPa A^3)
    @property
    def solubility_factor(self):
        if self.__solubility_factor_set == 0:
                self.__solubility_factor = self.average_shear_modulus*self.average_effective_lattice_constant * np.average(np.power(self.difference_effective_lattice_constant,2))/self.Z_value            
        return self.__solubility_factor
    @solubility_factor.setter
    def solubility_factor(self,value):
        self.__solubility_factor = value
        self.__solubility_factor_set = 1
        
    # 计算各组元比例分数
    @property
    def component_ratio(self):
        if self.component_proportion.size != 0:
        # 根据各个成分的比例计算摩尔分数
            self.__component_ratio = self.component_proportion/np.sum(self.component_proportion)
        return self.__component_ratio
   
    #计算构型熵(单位J/mol/K)
    @property
    def config_entropy(self):
        self.__config_entropy= float(-1.0 * R * np.sum(self.component_ratio * np.log(self.component_ratio)))     
        return self.__config_entropy
    #计算构型熵(单位R)
    @property
    def config_entropy_R(self):
        self.__config_entropy_R = float(-1.0 * np.sum(self.component_ratio * np.log(self.component_ratio)))
        return self.__config_entropy_R
    
    #计算熵变项(单位: J/mol)
    @property
    def entropy_term(self):
        self.__entropy_term = float(self.config_entropy * self.temperature)
        return self.__entropy_term
    #计算焓变项(单位: J/mol)
    @property
    def mixing_enthalpy(self):
        if self.component_number == 2:
            self.__mixing_enthalpy = float(self.M_value*N_A*self.solubility_factor*np.prod(self.component_ratio))
        elif (self.component_number > 2) and len(np.unique(self.component_ratio)==1):
            self.__mixing_enthalpy = float(self.M_value * N_A * self.solubility_factor*np.sum([(1-1/i)*(1/i)*np.power(i/self.component_number,3.5) for i in range(self.component_number+1)[2:]]))
        else:
            self.__mixing_enthalpy = 0
        return self.__mixing_enthalpy
    
    
    #计算吉布斯自由能(单位: J/mol)
    @property
    def mixing_Gibbus(self):
        self.__mixing_Gibbus = self.mixing_enthalpy - self.entropy_term
        return self.__mixing_Gibbus
    
def read_inpuut():
    prompt = "--->>>"
    CurrentPath=os.getcwd()
    AllFiles = os.listdir(CurrentPath)
    DataFile = []
    FileNum = 0
    for f in AllFiles:
        if os.path.splitext(f)[-1] == ".yaml":
            FileNum += 1
            print(f)
    if FileNum == 0:
        print("请准备输入文件。格式:xxx.yaml")
        print("""
            #必填设置
            Sample_Name: Bi0.3Sb1.7Te #样品名称，无实际意义，仅为区分
            Component_Name: [GeTe,MnTe,SnTe]
            Component_Proportion:[0.5,0.25,0.25]
            Lattice_Constant: [1.08,1.06,1.09]
            Shear_Modulus: [40,60,80]
            Z_Value: 4    # 结构基元数
            """)
        try:
            sys.exit(0)
        except:
            print("die")
        finally:
            print("clean up.")
            print("请准备输入文件。格式:xxx.yaml")
    YamlFile = os.path.join(CurrentPath,str(input(prompt)))
    with open(YamlFile,"r",encoding="utf-8") as f:
        value = yaml.load(f,Loader=yaml.FullLoader)
    return (value,YamlFile)
def calculate():
    print("-"*36+"基本信息"+"-"*36)
    print(" "*4+"本程序由13skeleton编写,如有任何问题，请直接联系邮箱。(J.Pei@foxmail.com)")
    print("""    参考文献：
    1. R. Liu, H. Chen, K. Zhao, Y. Qin, B. Jiang, T. Zhang, G. Sha, X. Shi, C. Uher, 
       W. Zhang, L. Chen, Entropy as a Gene-Like Performance Indicator Promoting Ther
       -moelectric Materials. Adv. Mater. 29, 1702712 (2017).
    2. Q. Yang, P. Qiu, X. Shi, L. Chen, Application of Entropy Engineering in Thermo
       -electrics. Journal of Inorganic Materials. 36, 347 (2021).
    3. M. Guo, F. Zhang, Y. Miao, Y. Liu, J. Yu, F. Gao, Preparation and Electrical 
       Properties of High Entropy La(Co 0.2 Cr 0.2 Fe 0.2 Mn 0.2 Ni 0.2 )O 3 Perovs
       -kite Ceramics Powder. Journal of Inorganic Materials. 36, 431 (2021).
      """)
    print("-"*36+"读取输入文件"+"-"*36)
    parameter,YamlFile = read_inpuut()
    s = High_Entropy(parameter["Sample_Name"])
    print("样品名称",s.sample_name)
    
    #读取输入文件
    if "Component_Name" in parameter.keys():
        s.component_name = np.array(parameter["Component_Name"])
        print("样品名称",s.component_name)
    if "Component_Number" in parameter.keys():
        s.component_number = np.array(parameter["Component_Number"])
    if "Component_Proportion" in parameter.keys():
        if len(parameter["Component_Proportion"]) == s.component_number:
            s.component_proportion = np.array(parameter["Component_Proportion"],dtype=np.float)
            print("各组元比例",s.component_proportion)
        else:
            print("各组元比例错误")
            exit(0)
    if "Lattice_Type" in parameter.keys():
        s.lattice_type = parameter["Lattice_Type"]
        print("晶格类型",s.lattice_type)
    if "Lattice_Constant" in parameter.keys():
        s.lattice_constant = np.array(parameter["Lattice_Constant"],dtype=np.float)
        print("晶格类型",s.lattice_constant)
    if "Effective_Lattice_Constant" in parameter.keys():
        s.effective_lattice_constant = np.array(parameter["Effective_Lattice_Constant"],dtype=np.float)
        print("晶格类型",s.effective_lattice_constant)
    if "Shear_Modulus" in parameter.keys():
        s.shear_modulus = np.array(parameter["Shear_Modulus"],dtype=np.float)
        print("剪切模量",s.shear_modulus)
    if "Average_Shear_Modulus" in parameter.keys():
        s.average_shear_modulus = float(parameter["Average_Shear_Modulus"])
        print("平均剪切模量",s.average_shear_modulus)
    if "Z_Value" in parameter.keys():
        s.Z_value = float(parameter["Z_Value"])
        print("Z值",s.Z_value)
    if "M_Value" in parameter.keys():
        s.M_value= float(parameter["M_Value"])
        print("无量纲因子M值",s.M_value)
    if "Solubility_Factor" in parameter.keys():
        s.solubility_factor = float(parameter["Solubility_Factor"])
        print("溶解度因子",s.solubility_factor)
    if "Temperature" in parameter.keys() and (not isinstance(parameter["Temperature"],dict)):
        s.temperature = float(parameter["Temperature"])
        print("温度",s.temperature)
        results_temperature = [s.temperature,s.entropy_term,s.mixing_Gibbus]

    if "Temperature" in parameter.keys() and isinstance(parameter["Temperature"],dict):        
        start_temperature = float(parameter["Temperature"]["Start_Temperature"])
        end_temperature = float(parameter["Temperature"]["End_Temperature"])
        interval_temperature = float(parameter["Temperature"]["Interval_Temperature"])
        print("起始温度",start_temperature)
        print("温度间隔",interval_temperature)
        print("终止温度",end_temperature)
        results_temperature = []
        for i in np.arange(start_temperature,end_temperature+interval_temperature,interval_temperature):
            s.temperature = i
            list_for_temperature = [s.temperature,s.entropy_term,s.mixing_Gibbus]
            results_temperature.append(list_for_temperature)
    print(" ")
    print(" ")
    print(" ")
    print("-"*36+"#输出结果"+"-"*36)
    if "Component_Number" not in parameter.keys():
        print("组分数量:",s.component_number)
    
    print("组分摩尔分数:",s.component_ratio)
    if "Lattice_Type" not in parameter.keys():
        print("基体的点阵类型:",s.lattice_type)
    if "Solubility_Factor" not in parameter.keys():                                         
        if "Effective_Lattice_Constant" not in parameter.keys():                                
            if "Lattice_Constant" in parameter.keys():
                print("有效晶胞常数:",s.effective_lattice_constant)
                print("晶胞常数之差:",s.difference_effective_lattice_constant)
                print("平均有效晶胞常数:",s.average_effective_lattice_constant)
        elif "Effective_Lattice_Constant" in parameter.keys():                                           
            print("晶胞常数之差:",s.difference_effective_lattice_constant)   
            print("平均有效晶胞常数:",s.average_effective_lattice_constant) 
        if "Average_Shear_Modulus" not in parameter.keys():
            if "Shear_Modulus" in parameter.keys():
                print("平均剪切模量:",s.average_shear_modulus)              

    if "Z_Value" not in parameter.keys():
        print("单胞中基元数目:",s.Z_value)

    if "Solubility_Factor" not in parameter.keys():
        print("溶解度因子:",s.solubility_factor)

    if "M_Value" not in parameter.keys():
        print("无量纲因子M:",s.M_value)

    print("构型熵",s.config_entropy)
    print("构型熵/R",s.config_entropy_R)
    print("混合焓",s.mixing_enthalpy)
    
    OutFile=os.getcwd()+os.sep+os.path.splitext(os.path.basename(YamlFile))[0]+"-out"+".csv"
    with open(OutFile,"w",newline="") as csvfile:
        myinput = csv.writer(csvfile)
        myinput.writerow(["#input parameter"])
        myinput.writerow(["sample name:",s.sample_name])
        #输入参数写入
        if "Component_Name" in parameter.keys():
            myinput.writerow(["样品名称",s.component_name])
        if "Component_Proportion" in parameter.keys():
            myinput.writerow(["各组元比例",s.component_proportion])
        if "Lattice_Type" in parameter.keys():
            myinput.writerow(["晶格类型",s.lattice_type])
        if "Lattice_Constant" in parameter.keys():
            myinput.writerow(["晶胞常数",s.lattice_constant])
        if "Effective_Lattice_Constant" in parameter.keys():
            myinput.writerow(["有效晶胞常数",s.effective_lattice_constant])
        if "Shear_Modulus" in parameter.keys():
            myinput.writerow(["剪切模量",s.shear_modulus])
        if "Average_Shear_Modulus" in parameter.keys():
            myinput.writerow(["平均剪切模量",s.average_shear_modulus])
        if "Z_Value" in parameter.keys():
            myinput.writerow(["Z值",s.Z_value])
        if "M_Value" in parameter.keys():
            myinput.writerow(["无量纲因子M值",s.M_value])
        if "Solubility_Factor" in parameter.keys():
            myinput.writerow(["溶解度因子",s.solubility_factor])
        if "Temperature" in parameter.keys() and (not isinstance(parameter["Temperature"],dict)):
            myinput.writerow(["温度",s.temperature])

        if "Temperature" in parameter.keys() and (isinstance(parameter["Temperature"],dict)):        
            myinput.writerow(["起始温度",start_temperature])
            myinput.writerow(["温度间隔",interval_temperature])
            myinput.writerow(["终止温度",end_temperature])
        myinput.writerow([" "," "])
        myinput.writerow([" "," "])
        
        myinput.writerow(["#output results"])
        if "Component_Number" not in parameter.keys():
            myinput.writerow(["组分数量:",s.component_number])

        myinput.writerow(["组分摩尔分数:",s.component_ratio])
        if "Lattice_Type" not in parameter.keys():
            myinput.writerow(["基体的点阵类型:",s.lattice_type])
        if "Solubility_Factor" not in parameter.keys():                                         
            if "Effective_Lattice_Constant" not in parameter.keys():                                
                if "Lattice_Constant" in parameter.keys():
                    myinput.writerow(["有效晶胞常数:",s.effective_lattice_constant])
                    myinput.writerow(["晶胞常数之差:",s.difference_effective_lattice_constant])
                    myinput.writerow(["平均有效晶胞常数:",s.average_effective_lattice_constant])
            elif "Effective_Lattice_Constant" in parameter.keys():                                           
                myinput.writerow(["晶胞常数之差:",s.difference_effective_lattice_constant])   
                myinput.writerow(["平均有效晶胞常数:",s.average_effective_lattice_constant]) 
            if "Average_Shear_Modulus" not in parameter.keys():
                if "Shear_Modulus" in parameter.keys():
                    myinput.writerow(["平均剪切模量:",s.average_shear_modulus])              

        if "Z_Value" not in parameter.keys():
            myinput.writerow(["单胞中基元数目:",s.Z_value])

        if "Solubility_Factor" not in parameter.keys():
            myinput.writerow(["溶解度因子:",s.solubility_factor])

        if "M_Value" not in parameter.keys():
            myinput.writerow(["无量纲因子M:",s.M_value])

        myinput.writerow(["构型熵",s.config_entropy])
        myinput.writerow(["构型熵/R",s.config_entropy_R])
        myinput.writerow(["混合焓",s.mixing_enthalpy])        
        
        myinput.writerow(["温度","熵变项","吉布斯自由能"])
        if "Temperature" in parameter.keys() and (not isinstance(parameter["Temperature"],dict)):
            myinput.writerow(results_temperature)
        if "Temperature" in parameter.keys() and (isinstance(parameter["Temperature"],dict)):
            myinput.writerows(results_temperature)
    print("...")
    print("...")
    print("...")
    
    print("计算完成，输出结果请查看xxx.csv文件")
if __name__ == "__main__":
    calculate()
    a=input("按任意键退出")