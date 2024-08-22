'''
Author: Jarvis-Chiang 497694894@qq.com
Date: 2024-08-10 16:19:57
LastEditors: Jarvis-Chiang 497694894@qq.com
LastEditTime: 2024-08-11 14:54:09
FilePath: /fenicsprojects/让我测试一下/save_data.py
Description: 从程序中获取处理后的数据
'''
import numpy as np
from dolfin import *

def save_vector_field_to_txt(mesh:mesh, theta_function:Function, filename="output/theta_vector_field.txt"):
    """
    将生成的向量场和网格数据保存成文本文件 (V-F 表)
    参数：
    mesh : dolfin.Mesh 
            有限元网格
    theta_function : dolfin.Function
            包含向量场信息的函数
    filename : str, optional 
            保存文件的名称，默认为 'theta_vector_field.txt'
    """

    

    # 获取网格节点坐标 (V 表)
    node_coordinates = mesh.coordinates()

    # 获取单元 (面) 的顶点索引 (F 表)
    faces = np.array([cell.entities(0) for cell in cells(mesh)])

    # 获取网格中心点坐标
    cell_midpoints = np.array([cell.midpoint().array() for cell in cells(mesh)])

    # 获取优化后的theta值
    theta_values = theta_function.vector().get_local()

    # 将角度转换为弧度
    theta_rad = np.deg2rad(theta_values)

    # 计算 x 和 y 方向的向量分量
    x_values = np.cos(theta_rad)
    y_values = np.sin(theta_rad)

    # 将坐标和向量分量组合成数组
    vector_field = np.column_stack((x_values, y_values))

#     # 保存到文本文件
#     np.savetxt(filename, vector_field, header="x y z u v", comments='', fmt='%.6f')
#     print(f"Vector field saved to {filename}")

    try:
        with open(filename, 'w') as f:
            
            # 写入 V 表 (顶点坐标)
            f.write("# ---------- Vertices (V table) ---------- \n")
            np.savetxt(f, node_coordinates, fmt='%.6f')      
            # 写入 F 表 (单元顶点索引)
            f.write("# ---------- Faces (F table) ---------- \n")
            np.savetxt(f, faces, fmt='%d')
            # 写入 U 表 (向量场)
            f.write("# ---------- Vector field (U table) ---------- \n")
            np.savetxt(f, vector_field, fmt='%.6f')

            print(f"File '{filename}' created successfully.")
            
    except IOError:
        print("An error occurred while trying to read the file.")
        return

