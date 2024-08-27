import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

class ShapeLoader:
    def __init__(self, file_path):
        # 初始化加载器，设置文件路径和数据结构
        self.file_path = file_path
        self.nodes = {}  # 存储节点信息
        self.elements = []  # 存储单元信息

    def load(self):
        # 从文件中加载数据
        with open(self.file_path, 'r') as file:
            current_section = None
            for line in file:
                line = line.strip()
                # 根据关键字确定当前处理的部分
                # 根据行首的标记来确定当前处理的部分
                if line.startswith('*NODE'):
                    current_section = 'node'  # 如果行首是'*NODE'，则当前处理的部分是节点
                elif line.startswith('*ELEMENT'):
                    current_section = 'element'  # 如果行首是'*ELEMENT'，则当前处理的部分是单元
                elif line.startswith('*'):
                    current_section = None  # 如果行首是'*'但不是'*NODE'或'*ELEMENT'，则当前处理的部分为None
                # 根据当前处理的部分来解析行
                elif current_section == 'node':
                    self._parse_node(line)  # 如果当前处理的部分是节点，则解析节点数据
                elif current_section == 'element':
                    self._parse_element(line)  # 如果当前处理的部分是单元，则解析单元数据

    def _parse_node(self, line):
        # 解析节点数据
        parts = line.split(',')
        node_id = int(parts[0])
        x, y = map(float, parts[1:3])  # 只取x和y坐标
        self.nodes[node_id] = (x, y)

    def _parse_element(self, line):
        # 解析单元数据
        parts = line.split(',')
        element_id = int(parts[0])
        node_ids = list(map(int, parts[1:]))
        self.elements.append((element_id, node_ids))

    def get_nodes(self):
        # 返回所有节点数据
        return self.nodes

    def get_elements(self):
        # 返回所有单元数据
        return self.elements

    def save_to_dolfin_xml(self, output_file):
        # 将数据保存为DOLFIN XML格式
        root = ET.Element("dolfin")
        
        mesh = ET.SubElement(root, "mesh", celltype="triangle", dim="2")  # 改为2D
        
        # 添加节点
        vertices = ET.SubElement(mesh, "vertices", size=str(len(self.nodes)))
        for node_id, (x, y) in self.nodes.items():
            ET.SubElement(vertices, "vertex", index=str(node_id-1), x=str(x), y=str(y))  # 只使用x和y
        
        # 添加单元
        cells = ET.SubElement(mesh, "cells", size=str(len(self.elements)))
        for i, (element_id, node_ids) in enumerate(self.elements):
            ET.SubElement(cells, "triangle", index=str(i), v0=str(node_ids[0]-1), v1=str(node_ids[1]-1), v2=str(node_ids[2]-1))
        
        # 美化XML输出
        xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
        
        # 写入文件
        with open(output_file, "w") as f:
            f.write(xml_str)
    
if __name__ == "__main__":
    # 主程序入口
    loader = ShapeLoader("input/topo_shape.inp")
    loader.load()
    # print(loader.get_nodes())
    # print(loader.get_elements())
    loader.save_to_dolfin_xml("input/topo_shape.xml")