'''
Author: Jarvis-Chiang 497694894@qq.com
Date: 2024-08-10 18:50:32
LastEditors: Jarvis-Chiang 497694894@qq.com
LastEditTime: 2024-08-10 18:53:07
FilePath: /fenicsprojects/让我测试一下/show_stream_line.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import pyvista as pv 
import numpy as np

def plot_vector_field(grid, **options):
    """
    使用 PyVista 绘制向量场，包括箭头和流场。

    参数:
    grid: pv.PolyData
        包含点和向量数据的PolyData对象。
    options: dict, optional
        可选参数，包括是否显示title、scalar_bar、label。
    """
    streamlines = grid.streamlines('vectors', n_points=200, max_time=2.0)
    plotter.add_mesh(streamlines, color="red")

    show_grid = options.get("show_grid", True)
    show_title = options.get("show_title", True)
    show_scalar_bar = options.get("show_scalar_bar", False)
    show_streamlines = options.get("show_streamlines", True)

    # 使用 glyph 方法创建箭头表示向量
    arrows = grid.glyph(orient="vectors", scale=False, factor=0.1)
        # 使用 PyVista 进行可视化
    plotter = pv.Plotter()

    if show_grid:
        plotter.show_bounds(grid='back', location='outer', all_edges=True)

    plotter.add_mesh(arrows, color="blue", show_scalar_bar=show_scalar_bar)

    if show_title:
        plotter.add_text("Vector Field Visualization", position="upper_edge", color="black", font_size=10)



    plotter.view_xy()
    plotter.show()