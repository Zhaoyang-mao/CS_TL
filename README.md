<p align="center">
  <br><strong>一种基于迁移学习的感知矩阵优化方法</strong> 
  <br>毛朝阳，李岚，魏伟
</p>

## 描述
#### 题目:
《一种基于迁移学习的感知矩阵优化方法》；《An optimization method of sensing matrix based on transfer learning》
#### 摘要
压缩感知中提高信号重构精度的关键问题是设计有效的感知矩阵。本文提出基于迁移学习的感知矩阵优化方法。首先通过迁移学习更新稀疏表示系数，将固定稀疏基转换为自适应的稀疏基。其次根据稀基与测量矩阵的乘积构造一个Gram矩阵，通过特征分解最小化其非对角线元素，减小Gram矩阵的全局相干性，从而实现原始信号的精确重建 。

## 文件结构
- **data**  包含测试数据。
- **Exmple** 测试实例，可直接运行。
- **OptMatrixMethod** 包含四种测量矩阵的优化算法
- **RecoverAlgorithm** 压缩感知的恢复算法，用的是**OMP** 算法。
- **Src** 包含一些其它辅助文件
