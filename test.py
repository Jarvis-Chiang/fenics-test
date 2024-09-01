'''
Author: Jarvis-Chiang 497694894@qq.com
Date: 2024-07-10 10:38:59
LastEditors: Jarvis-Chiang 497694894@qq.com
LastEditTime: 2024-08-02 11:21:49
FilePath: /fenicsprojects/让我测试一下/test.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%A9
'''
import vtk

def main():
    # 创建一个 Sphere 对象
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(5.0)
    sphere.SetPhiResolution(20)
    sphere.SetThetaResolution(20)

    # 创建一个 Mapper 对象
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())

    # 创建一个 Actor 对象
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    # 添加材质
    material = vtk.vtkProperty()
    material.SetColor(1.0, 0.5, 0.31)
    material.SetSpecular(0.5)
    material.SetSpecularPower(30)
    material.SetOpacity(0.5)
    actor.SetProperty(material)

    # 显示实心网格而不是线框
    actor.GetProperty().EdgeVisibilityOff()
    actor.GetProperty().SetRepresentationToSurface()

    # 显示网格实体
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(1.0, 1.0, 1.0)

    # 创建坐标系
    axes = vtk.vtkAxesActor()
    axes.SetTotalLength(1.5, 1.5, 1.5)
    axes.AxisLabelsOff()
    
    # 保持坐标系的原点不动，不进行平移
    transform = vtk.vtkTransform()
    transform.Translate(-7.5, -7.5, 0)
    axes.SetUserTransform(transform)

    # 创建一个 Renderer 对象
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor(axes)
    renderer.SetBackground(0.9, 0.9, 0.9)

    # 创建一个 RenderWindow 对象
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # 创建一个 RenderWindowInteractor 对象
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # 开始渲染并交互
    renderWindow.Render()
    renderWindowInteractor.Start()

if __name__ == "__main__":
    main()
