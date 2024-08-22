'''
Author: Jarvis-Chiang 497694894@qq.com
Date: 2024-08-10 16:58:31
LastEditors: Jarvis-Chiang 497694894@qq.com
LastEditTime: 2024-08-16 11:16:19
FilePath: ./show_vector_with_pyvista.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import pyvista as pv
import numpy as np


class ShowWithPyvista:
    def __init__(self, filename, **options):
        """
        初始化类并读取向量场数据。

        参数:
        filename : str
            包含向量场数据的文本文件的名称。
        """
        self.vertices = []
        self.faces = []
        self.midpoints = []        
        self.vector_field_midpints = []
        self.vector_field = []
        self.grid = pv.UnstructuredGrid()
        self.plotter = pv.Plotter(**options)
        self.stream_start_points = []
        
        self._set_plotter(view_xy=True)
        self._load_data(filename)
        self._calculate_midpoints()
        self._interpolate_vector_field()
        self._creat_grid()


    def _set_plotter(self, **options):
        # 获取设置参数
        plotter_title = options.get("title", "Plotter Title")   # 设置图形title
        view_xy = options.get("view_xy", False)     # xy视角

        self.plotter.add_text(plotter_title, position='upper_edge')
        if view_xy:
            self.plotter.view_xy()

    def _load_data(self, filename):
        """
        从文件中加载向量场数据。

        参数:
        filename : str
            包含向量场数据的文本文件的名称。
        """
        try:
            with open(filename, 'r') as file:
                content = file.readlines()
        except FileNotFoundError:
            print("File not found. Please input a valid file name!")
            return
        except IOError:
            print("An error occurred while trying to read the file.")
            return

        current_section = None
        for line in content:
            if line.startswith("#"):
                if "Vertices" in line:
                    current_section = "vertices"
                elif "Faces" in line:
                    current_section = "faces"
                elif "Vector field" in line:
                    current_section = "vector_field"
                continue
            if line.strip():
                data = list(map(float, line.strip().split()))
                if current_section == "vertices":
                    self.vertices.append(data)
                elif current_section == "faces":
                    self.faces.append([int(x) for x in data])
                elif current_section == "vector_field":
                    self.vector_field_midpints.append(data)

        # 填充 self.vertices 变量
        self.vertices = np.array(self.vertices)
        if self.vertices.shape[1] == 2:
            self.vertices = np.hstack((self.vertices, np.zeros_like(self.vertices[:, 0]).reshape(-1, 1)))

        # 填充 self.faces 变量
        self.faces = np.array(self.faces, dtype=int)

        # 填充 self.vector_field 变量
        self.vector_field_midpints = np.array(self.vector_field_midpints)
        if self.vector_field_midpints.shape[1] == 2:
            self.vector_field_midpints = np.hstack((self.vector_field_midpints, np.zeros_like(self.vector_field_midpints[:, 0]).reshape(-1, 1)))

        # 填充 self.midpoints 变量
        self._calculate_midpoints()

        # 填充 self.grid 变量
        self._creat_grid()

    def __str__(self):
        # 计算和返回需要的统计信息
        faces_count = len(self.faces)
        vertex_count = len(self.vertices)
        vector_field_count = len(self.vector_field)
        return (f"Grid Count: {faces_count}\n"
                f"Vertex Count: {vertex_count}\n"
                f"Vector Field Count: {vector_field_count}")

    def _calculate_midpoints(self):
        """
        计算每个面的中心点，并将其存储在 self.midpoints 中。
        """
        self.midpoints = []
        for face in self.faces:
            face_vertices = self.vertices[face]
            midpoint = np.mean(face_vertices, axis=0)
            self.midpoints.append(midpoint)
        self.midpoints = np.array(self.midpoints)

    def _creat_grid(self):
        cells = grid_faces = np.insert(self.faces, 0, 3, axis=1) 
        celltypes = np.full(len(cells), pv.CellType.TRIANGLE)
        points = self.vertices
        self.grid = pv.UnstructuredGrid(cells, celltypes, points)

    def _interpolate_vector_field(self, power=2):
        """
        使用逆距离加权插值（IDW）将中心点的向量场插值到网格顶点上。

        参数:
        power : float
            距离的权重指数，默认值为2。
        """
        self.vector_field = np.zeros_like(self.vertices)

        for i, vertex in enumerate(self.vertices):
            # 计算从该顶点到所有中心点的距离
            distances = np.linalg.norm(self.midpoints - vertex, axis=1)
            
            # 防止距离为零，添加一个小值
            distances[distances == 0] = 1e-10
            
            # 计算权重
            weights = 1.0 / distances**power
            
            # 计算加权平均
            weighted_vectors = np.sum(weights[:, np.newaxis] * self.vector_field_midpints, axis=0)
            self.vector_field[i] = weighted_vectors / np.sum(weights)        

    def calculate_vector_at_point(p, vertices, vector_field, power=2):
        """
        计算网格内任一点的向量大小和方向。

        参数:
         @p : array-like
            待插值点的坐标。
         @vertices : array-like
            网格顶点的坐标。
         @vector_field : array-like
            网格顶点处的向量场。
         @power : float
            距离的权重指数，默认值为2。

        返回:
         @interpolated_vector : array
            插值点处的向量。
         @magnitude : float
            插值点处的向量大小。
         @direction : array
            插值点处的向量方向（归一化向量）。
        """
        p = np.array(p)
        distances = np.linalg.norm(vertices - p, axis=1)
        distances[distances == 0] = 1e-10  # 防止距离为零，添加一个小值

        # 计算权重
        weights = 1.0 / distances**power

        # 计算加权平均
        interpolated_vector = np.sum(weights[:, np.newaxis] * vector_field, axis=0) / np.sum(weights)

        # 计算向量大小和方向
        magnitude = np.linalg.norm(interpolated_vector)
        direction = interpolated_vector / magnitude if magnitude != 0 else np.zeros(3)

        return interpolated_vector, magnitude, direction

    def add_vector_field_midpoint(self, **options):
        """
        使用 PyVista 显示向量场数据。

        参数:
        **options : dict
            动态选项参数，包括是否显示网格、标题、标量条、标签等。
            可选参数:
            - show_grid: bool, 是否显示网格，默认值为 True。
            - show_title: bool, 是否显示标题，默认值为 False。
            - show_scalar_bar: bool, 是否显示标量条，默认值为 False。
            - show_label: bool, 是否显示标签，默认值为 False。
        """
        show_mesh = options.get('show_mesh', True)
        title = options.get('title', "vector field")
        show_scalar_bar = options.get('show_scalar_bar', False)
        show_label = options.get('show_label', False)
        scalar = options.get('scalar', 0.5)

        # 确保 midpoints 被计算
        if not hasattr(self, 'midpoints') or len(self.midpoints) == 0:
            self._calculate_midpoints()

        # 将向量场绘制在 midpoints 上
        # 检测维度是否相同
        if len(self.midpoints) != len(self.vector_field_midpints):
            raise ValueError("The length of midpoints and vector_field must be the same.")
        
        self.plotter.add_arrows(self.midpoints, self.vector_field_midpints[:len(self.midpoints)], mag=scalar)
        return
    
    def add_vector_field(self, **options):

        # 获取设置参数
        mag = options.get("mag", 0.3)
        scalars = options.get("scalars", 1.0)
        
        # # 检测 vector_field 维度是否和网格顶点个数相同
        # if len(self.vector_field) != len(self.vertices):
        #     raise ValueError("The number of vector field and vertices must be the same.")

        # # 如果 scalars 数组长度不匹配，进行扩展
        # scalars = np.full(len(self.vertices), scalars)
        
        # 添加箭头
        self.plotter.add_arrows(self.vertices, self.vector_field, mag=mag)
        return

    def add_start_points(self, start_points):
        """
        给流线图添加流线起点

        参数：
        - start_point : array like
            流线图的起点坐标
        """
        self.stream_start_points = start_points
        return 

    def add_streamlines(self, start_points=[], **options):
        """
        使用 PyVista 显示流线图。

        参数:
        start_points : array-like
            流线图的起始点坐标。
        **options : dict
            动态选项参数，包括流线图的其他设置。
            可选参数:
            - show_grid: bool, 是否显示网格，默认值为 True。
            - show_title: bool, 是否显示标题，默认值为 False。
            - show_scalar_bar: bool, 是否显示标量条，默认值为 False。
        """

        # 检测是否满足前置条件
        if not hasattr(self, 'grid'):
            raise AttributeError("Grid is not defined in the object.")
        
        self.grid['vectors'] = self.vector_field
        
        if len(start_points) == 0:
            raise ValueError("No start points provided for streamlines.")
        
        # 添加起始点
        self.add_start_points(start_points)
        # 将起始点转换为 PyVista 的 PolyData 对象
        source = pv.PolyData(self.stream_start_points)
        
        streamlines = self.grid.streamlines_from_source(
            source=source,
            vectors='vectors',  # 指定向量场的名称
            integration_direction='both',  # 流线的方向，'both'表示在正反两个方向上生成流线
            max_time=50.0,  # 控制流线的长度，值越大，流线越长
            max_steps=20000,  # 流线生成的最大步数，步数越多，流线越长
            terminal_speed=1e-12,  # 当速度小于该值时，终止流线生成
            initial_step_length=0.01,  # 流线生成的初始步长
            max_step_length=0.1,  # 流线生成的最大步长，控制流线延展的最大步伐
            min_step_length=0.1,  # 流线生成的最小步长，控制流线延展的最小步伐
            step_unit='cl',  # 步长的单位，'cl' 表示以单元长度为单位
            progress_bar=True,
        )
        
        # 检测是否生成了流线
        if streamlines.n_points == 0:
            raise RuntimeError("No streamlines were generated.")
        
        self.plotter.add_mesh(streamlines.tube(radius=0.01), color="black")
        
        

        if options.get('show_scalar_bar', False):
            self.plotter.add_scalar_bar()

        return




        # if len(start_points) != 0:
        #     self.stream_start_points = start_points

        # if len(self.stream_start_points) == 0 and len(start_points) == 0:
        #     raise ValueError("Must have start point.")

        # show_grid = options.get('show_grid', True)
        # title = options.get('title', "Streamlines")
        # show_scalar_bar = options.get('show_scalar_bar', False)

        # # 确保 midpoints 被计算
        # if not hasattr(self, 'midpoints') or len(self.midpoints) == 0:
        #     self._calculate_midpoints()

        # # 使用 PyVista 进行可视化
        # grid_faces = np.insert(self.faces, 0, 3, axis=1) 
        # grid = pv.PolyData(self.vertices, grid_faces)

        # # 创建一个 PolyData 对象用于流线图
        # vector_field = np.array(self.vector_field_midpints)
        # grid.point_data['vectors'] = vector_field

        # self.plotter.add_mesh(grid, show_edges=show_grid, opacity=0.1)

        # # 添加流线图
        # streamlines = grid.streamlines(start_points=start_points, vectors='vectors')
        # self.plotter.add_mesh(streamlines, color='blue')

        # if show_grid:
        #     self.plotter.add_mesh(grid, show_edges=show_grid, opacity=0.2)

        # if show_scalar_bar:
        #     self.plotter.add_scalar_bar("Magnitude")

        # self.plotter.add_title(title)

    def add_grid(self, **options):
        self.plotter.add_mesh(self.grid)


    def show(self):
        self.plotter.show()

if __name__ == "__main__":
    sp = ShowWithPyvista("/home/wsl-20/fenicsprojects/让我测试一下/output/theta_vector(adam).txt")
    print("pause!!")