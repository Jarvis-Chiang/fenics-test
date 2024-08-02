'''
Author: Jarvis-Chiang 497694894@qq.com
Date: 2024-07-10 10:38:59
LastEditors: Jarvis-Chiang 497694894@qq.com
LastEditTime: 2024-08-02 11:21:49
FilePath: /fenicsprojects/让我测试一下/test.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import vtk

# 创建一个 Sphere 对象
sphere = vtk.vtkSphereSource()
sphere.SetRadius(5.0)
sphere.SetPhiResolution(50)
sphere.SetThetaResolution(50)

# 创建一个 Mapper 对象
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(sphere.GetOutputPort())

# 创建一个 Actor 对象
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# 创建一个 Renderer 对象
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(0.1, 0.2, 0.4)  # 背景颜色

# 创建一个 RenderWindow 对象
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# 创建一个 RenderWindowInteractor 对象
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# 开始渲染并交互
renderWindow.Render()
renderWindowInteractor.Start()