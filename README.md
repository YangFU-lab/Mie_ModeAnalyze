[README_Mie_ModeAnalysis.md](https://github.com/user-attachments/files/23713372/README_Mie_ModeAnalysis.md)
# Mie散射模式分解分析

基于Mie散射理论，分解并分析球形粒子散射中各个电偶极和磁偶极模式的贡献。

## 作者信息

- **姓名**: Yang FU
- **研究网站**: [https://yfu-research.com/](https://yfu-research.com/)

## 引用

如果您在研究中使用了本代码，请引用以下文献，十分感谢：

**1. Fu, Y., An, Y., Xu, Y., Dai, J. G., & Lei, D. (2022). Polymer coating with gradient‐dispersed dielectric nanoparticles for enhanced daytime radiative cooling. *EcoMat*, 4(2), e12169.**  
https://doi.org/10.1002/eom2.12169

**2. Abou-Hamdan, L., Coudrat, L., Bidault, S., Krachmalnicoff, V., Haïdar, R., Bouchon, P., & De Wilde, Y. (2022). Transition from phononic to geometrical mie modes measured in single subwavelength polar dielectric spheres. *ACS photonics*, 9(7), 2295-2303.**  
https://doi.org/10.1021/acsphotonics.2c00273

## 功能

- 计算总散射截面和吸收截面
- 分解各Mie模式（an和bn）的贡献
- 可视化各模式的散射和吸收截面
- 支持自定义波长范围和粒子参数

## 依赖项

### 必需函数
- `Mie.m` - Mie散射计算函数
- `Mie_ab.m` - Mie系数计算函数

### 材料光学常数函数
- `sio2_n(lamda)` - SiO₂折射率
- `sio2_k(lamda)` - SiO₂消光系数

可根据需要修改为其他材料函数。

## 使用方法

1. 确保依赖函数在MATLAB路径中
2. 运行 `Mie_ModeAnalysis.m`
3. 根据需要修改参数（见下方）

## 主要参数

### 默认设置
- **波长范围**: 2.5 - 30 μm，步长 0.01 μm
- **粒子半径**: 25 nm
- **分析模式数**: 5（an₁-an₅, bn₁-bn₅）
- **介质**: 空气（n=1, k=0）

### 修改参数示例

```matlab
% 修改波长范围
lamda = (0.4:0.01:0.7)' * 10^-6;  % 可见光范围

% 修改粒子半径
r_particle = 100 * 1e-9;  % 100 nm

% 修改分析模式数
N = 3;  % 只分析前3个模式

% 修改材料
n_particle = TiO2_n(lamda);
k_particle = TiO2_k(lamda);
```

## 输出结果

### 计算结果
- `Cscat` - 总散射截面 [波长 × 1]
- `Cabs` - 总吸收截面 [波长 × 1]
- `Cscat_an` - 各an模式散射截面 [波长 × 模式数]
- `Cscat_bn` - 各bn模式散射截面 [波长 × 模式数]
- `Cabs_an` - 各an模式吸收截面 [波长 × 模式数]
- `Cabs_bn` - 各bn模式吸收截面 [波长 × 模式数]

### 图形输出
1. 介质光学常数图
2. 粒子光学常数图
3. 粒子介电常数图
4. 散射截面模式分解图
5. 吸收截面模式分解图

## Mie模式说明

- **an模式** (n=1,2,3...): 电偶极、电四极、电八极等
- **bn模式** (n=1,2,3...): 磁偶极、磁四极、磁八极等

总散射/吸收 = Σ(各模式贡献)

## 应用场景

- 分析不同尺寸粒子的散射特性
- 研究各模式对总散射的贡献
- 优化粒子尺寸以获得特定模式共振
- 理解Mie散射的物理机制

## 注意事项

- 计算时间取决于波长数量和模式数
- 小粒子主要贡献来自低阶模式（n=1,2）
- 大粒子需要更多模式才能准确计算

## 参考文献

- Bohren, C. F., & Huffman, D. R. (1983). *Absorption and Scattering of Light by Small Particles*. Wiley.
