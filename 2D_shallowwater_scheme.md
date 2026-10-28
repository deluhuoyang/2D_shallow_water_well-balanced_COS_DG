# 二维浅水波方程

## 1. well-balanced格式

二维浅水波方程：
$$
\begin{cases}
h_t + (hu)_x + (hv)_y = 0, \\
(hu)_t + \left(hu^2+\frac{1}{2}gh^2\right)_x + (huv)_y = -ghb_x, \\
(hv)_t + (huv)_x + \left(hv^2+\frac{1}{2}gh^2\right)_y = -ghb_y,
\end{cases}
$$

其中 \(h\) 为水深，\((u,v)\) 为流速，\(b(x,y)\) 为底部高程，\(g\) 为重力加速度。

为便于表述，我们将该方程记作
$$
U_t + \boldsymbol{f}(U)_x + \boldsymbol{g}(U)_y = s(h,b),
$$
其中 $U = (h,hu,hv)^T$，$\boldsymbol{f}(U),\boldsymbol{g}(U)$ 为通量。

我们将计算区域离散为单元
$$
I_{i,j} = \left[x_{i-\frac12},x_{i+\frac12}\right] \times \left[y_{j-\frac12},y_{j+\frac12}\right].
$$

为简化起见，采用均匀网格，记 $\lambda_1=\Delta t/\Delta x,\ \lambda_2=\Delta t/\Delta y$。在每个单元 $I_{i,j}$ 内，解空间与检验空间选取次数至多为 $k$ 的二元多项式空间。数值近似解为分片多项式，仍记为 $U$；底坡 $b$ 投影到有限元空间后的结果也仍记作 $b$。

与一维的情形理由相同，二维DG方法表述为：
$$
\begin{aligned}
\int_{I_{ij}} \frac{U^{n+1}-U^n}{\Delta t} v \,dxdy
&-\int_{I_{ij}} f(U^n)\partial_x v \,dxdy
+\int_{y_{j-\frac12}}^{y_{j+\frac12}} \hat{f}_{i+\frac12} v\left(x_{i+\frac12}^-,y\right) dy \\
&-\int_{y_{j-\frac12}}^{y_{j+\frac12}} \hat{f}_{i-\frac12} v\left(x_{i-\frac12}^+,y\right) dy
-\int_{I_{ij}} g(U^n)\partial_y v \,dxdy \\
&+\int_{x_{i-\frac12}}^{x_{i+\frac12}} \hat{g}_{j+\frac12} v\left(x,y_{j+\frac12}^-\right) dx
-\int_{x_{i-\frac12}}^{x_{i+\frac12}} \hat{g}_{j-\frac12} v\left(x,y_{j-\frac12}^+\right) dx \\
&=\int_{I_{ij}} s(h^n,b) v \,dxdy,
\end{aligned}
$$

其中 $\hat{f}_{i+\frac12}(y)=F\left(U\left(x_{i+\frac12}^-,y,t\right),U\left(x_{i+\frac12}^+,y,t\right)\right)$，采用Lax-Friedrichs通量：

$F(a_1,a_2)=\frac12\big(f(a_1)+f(a_2)-\alpha_1(a_2-a_1)\big),$

$\alpha_1=\max\left(|u|+\sqrt{gh}\right),$

且 $\hat{g}_{j+\frac12}(x)=G\left(U\left(x,y_{j+\frac12}^-\right),U\left(x,y_{j+\frac12}^+,t\right)\right)$，其中

$G(a_1,a_2)=\frac12\big(g(a_1)+g(a_2)-\alpha_2(a_2-a_1)\big),$

$\alpha_2=\max\left(|v|+\sqrt{gh}\right).$

将一维的保平衡RKDG格式推广到二维是直接可行的。设 $\hat{f}_{i+\frac12}^l,\hat{f}_{i-\frac12}^r,\hat{g}_{j+\frac12}^l,\hat{g}_{j-\frac12}^r$ ，则保平衡DG方法为
$$
\begin{aligned}
\int_{I_{ij}} \frac{U^{n+1}-U^n}{\Delta t} v \,dxdy
&-\int_{I_{ij}} f(U^n)\partial_x v \,dxdy
+\int_{y_{j-\frac12}}^{y_{j+\frac12}} \hat{f}_{i+\frac12}^l v\left(x_{i+\frac12}^-,y\right) dy \\
&-\int_{y_{j-\frac12}}^{y_{j+\frac12}} \hat{f}_{i-\frac12}^r v\left(x_{i-\frac12}^+,y\right) dy
-\int_{I_{ij}} g(U^n)\partial_y v \,dxdy \\
&+\int_{x_{i-\frac12}}^{x_{i+\frac12}} \hat{g}_{j+\frac12}^l v\left(x,y_{j+\frac12}^-\right) dx
-\int_{x_{i-\frac12}}^{x_{i+\frac12}} \hat{g}_{j-\frac12}^r v\left(x,y_{j-\frac12}^+\right) dx \\
&=\int_{I_{ij}} s(h^n,b) v \,dxdy,
\end{aligned}
$$

## 2. COS 过程

以二维矩形单元 $K_{i,j} = [x_{i-\frac12},x_{i+\frac12}] \times [y_{j-\frac12},y_{j+\frac12}]$ 为例，$\boldsymbol{x}=(x,y)$ 为空间坐标。相邻单元 $K_{i\pm1,j}$ 与 $K_{i,j\pm1}$ 分别和 $K_{i,j}$ 共用一条边。定义各坐标方向的空间步长：$h_{x,i}=x_{i+\frac12}-x_{i-\frac12}$，$h_{y,j}=y_{j+\frac12}-y_{j-\frac12}$。单元 $K_{i,j}$ 四条边上对应的局部最大波速定义为

$
\begin{aligned}
\beta_{i\pm\frac12,j}(\boldsymbol{u}_h^n) &= \max\left\{ r\left(\frac{\partial \boldsymbol{f}}{\partial \boldsymbol{u}}(\bar{\boldsymbol{u}}_{i\pm\frac12,j}^n)\right),\ r\left(\frac{\partial \boldsymbol{f}}{\partial \boldsymbol{u}}(\bar{\boldsymbol{u}}_{i,j}^n)\right) \right\},\\
\beta_{i,j\pm\frac12}(\boldsymbol{u}_h^n) &= \max\left\{ r\left(\frac{\partial \boldsymbol{g}}{\partial \boldsymbol{u}}(\bar{\boldsymbol{u}}_{i,j\pm\frac12}^n)\right),\ r\left(\frac{\partial \boldsymbol{g}}{\partial \boldsymbol{u}}(\bar{\boldsymbol{u}}_{i,j}^n)\right) \right\}.
\end{aligned}
$

于是系数 $\lambda_{K_{i,j}}(\boldsymbol{u}_h^n)$ 写为

$\lambda_{K_{i,j}}(\boldsymbol{u}_h^n)=\exp\left(-\Delta t\left(\frac{\beta_{i-\frac12,j}\sigma_{i-\frac12,j} + \beta_{i+\frac12,j}\sigma_{i+\frac12,j}}{h_{x,i}}
+\frac{\beta_{i,j-\frac12}\sigma_{i,j-\frac12} + \beta_{i,j+\frac12}\sigma_{i,j+\frac12}}{h_{y,j}}
\right)
\right)
$

其中
$
\sigma_{i\pm\frac12,j}(\boldsymbol{u}_h^n)=
\begin{cases}
C_k \dfrac{\|\boldsymbol{u}_{i\pm1,j}^n(\boldsymbol{x}) - \boldsymbol{u}_{i,j}^n(\boldsymbol{x})\|_{L_g^2(K_{i\pm1,j})}^2}{\Theta_{i\pm1,j}},
& \text{若 } \Theta_{i\pm1,j} \ge \varepsilon,\\[6pt]
0, & \text{其他情况},
\end{cases}
$

$
\sigma_{i,j\pm\frac12}(\boldsymbol{u}_h^n)=
\begin{cases}
C_k \dfrac{\|\boldsymbol{u}_{i,j\pm1}^n(\boldsymbol{x}) - \boldsymbol{u}_{i,j}^n(\boldsymbol{x})\|_{L_g^2(K_{i,j\pm1})}^2}{\Theta_{i,j\pm1}},
& \text{若 } \Theta_{i,j\pm1} \ge \varepsilon,\\[6pt]
0, & \text{其他情况},
\end{cases}
$

式中
$
\Theta_{i\pm1,j} = \|\boldsymbol{u}_{i\pm1,j}^n(\boldsymbol{x}) - \bar{\boldsymbol{u}}_{i,j}^n\|_{L_g^2(K_{i\pm1,j})}^2,\quad
\Theta_{i,j\pm1} = \|\boldsymbol{u}_{i,j\pm1}^n(\boldsymbol{x}) - \bar{\boldsymbol{u}}_{i,j}^n\|_{L_g^2(K_{i,j\pm1})}^2,
$

$
\varepsilon = 10^{-8}\max\big\{g(\bar{\boldsymbol{u}}_K^n),\,1\big\},
$

参数 $C_k$ 在本实验统一取为 0.2。

浅水波方程（SWE）的熵函数（能量）取为严格凸函数：
$$
g(\boldsymbol{u}) = \frac12\left( \frac{(hu)^2+(hv)^2}{h} + g h^2 \right)
$$

---

## 3. C-property 验证

### 3.1 初始条件

取矩形计算域：周期边界 $\Omega = [0,1]\times[0,1].$

底部：
$$
b(x,y)=0.8\exp\left(-50\left((x-0.5)^2+(y-0.5)^2\right)\right).
$$

初始水深为：
$$
h(x,y,0)=1-b(x,y),
$$

初始动量为：
$$
hu(x,y,0)=0,\qquad hv(x,y,0)=0.
$$

### 3.2 实验结果

**表2.2.1** 不同精度下稳态解的 $L^1$ 与 $L^\infty$ 误差（$mesh: 100\times100, CFL = 0.16, tend = 0.1 s$）

| Precision | $L^1$ error |  |  | $L^\infty$ error |  |  |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|  | $h$ | $hu$ | $hv$ | $h$ | $hu$ | $hv$ |
| Single | $1.56\text{e-07}$ | $9.23\text{e-07}$ | $9.55\text{e-07}$ | $7.75\text{e-07}$ | $5.52\text{e-06}$ | $6.50\text{e-06}$ |
| Double | $5.42\text{e-15}$ | $3.32\text{e-15}$ | $3.41\text{e-15}$ | $8.33\text{e-15}$ | $1.82\text{e-14}$ | $1.82\text{e-14}$ |

---

## 4. 收敛阶验证

### 4.1 初始条件
底地形：
$$
b(x,y) = \sin^2(\pi x)\sin^2(\pi y)
$$
或
$$
b(x,y) = \sin^2(\pi x) + \sin^2(\pi y)
$$

初始条件：
$$
h(x,y,0) = 5 + e^{\cos(2\pi x)} e^{\cos(2\pi y)}
$$
$$
(hu)(x,y,0) = \sin\big(\cos(2\pi x)\big),
$$
$$
(hv)(x,y,0) = \sin\big(\cos(2\pi y)\big)
$$

计算到 $t=0.05$。周期边界 $[0,1]^2$。

### 4.2 实验结果

**以128x128的网格为精确解。**（$CFL = 0.16, tend = 0.05 s$）
| Mesh | $h$ |  | $hu$ |  | $hv$ |  |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
|  | $L^1$ error | Order | $L^1$ error | Order | $L^1$ error | Order |
| $4\times4$ | $2.2001\text{e-01}$ | --- | $1.1485\text{e+00}$ | --- | $1.1485\text{e+00}$ | --- |
| $8\times8$ | $5.9420\text{e-02}$ | $1.89$ | $2.8607\text{e-01}$ | $2.01$ | $2.8607\text{e-01}$ | $2.01$ |
| $16\times16$ | $7.4332\text{e-03}$ | $3.00$ | $4.2694\text{e-02}$ | $2.74$ | $4.2694\text{e-02}$ | $2.74$ |
| $32\times32$ | $8.4310\text{e-04}$ | $3.14$ | $5.4796\text{e-03}$ | $2.96$ | $5.4796\text{e-03}$ | $2.96$ |
| $64\times64$ | $7.3461\text{e-05}$ | $3.52$ | $4.8299\text{e-04}$ | $3.50$ | $4.8299\text{e-04}$ | $3.50$ |


---

## 5. 二维振荡湖算例
### 5.1 初始条件
计算区域取矩形区域 \([-2,2]\times[-2,2]\)。抛物型底地形形式为：
$$
b(x,y)=h_0\frac{x^2+y^2}{a^2}
$$

其中常数 \(h_0,a\) 后续给定。对无摩擦源项的二维浅水波方程，解析解由论文给出：
$$
\begin{aligned}
h(x,y,t)&=\max\left(0,\frac{\sigma h_0}{a^2}\big(2x\cos(\omega t)+2y\sin(\omega t)-\sigma\big)+0.1-b(x,y)\right),\\
u(x,y,t)&=-\sigma\omega\sin(\omega t),\quad v(x,y,t)=\sigma\omega\cos(\omega t),
\end{aligned}
$$

该解以 \(T=2\pi/\omega\) 为周期，且 $\omega=\sqrt{2gh_0}/a$。

本算例固定参数：\(a=1,\,\sigma=0.5,\,h_0=0.1\)。初始条件由式 (6.13) 取 \(t=0\) 得到。由于流动不会抵达区域边界，边界条件可以任意选取，对数值解无影响。

我们采用 \(100\times100\) 的均匀网格，$CFL = 0.16$

### 5.2 实验结果
以下图片左上：$t = 0$，右上：$t = T/6$，左下：$t = T/3$，右下：$t = T/2$。

#### 5.2.1 3D 图
<div align="center">
<table>
  <tr>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_cell_100x100.png" width="400"/>
      <br/>
      $t=0$
    </td>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_T0.747541_cell_100x100.png" width="400"/>
      <br/>
      $t=T/6$
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_T1.495081_cell_100x100.png" width="400"/>
      <br/>
      $t=T/3$
    </td>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_T2.242622_cell_100x100.png" width="400"/>
      <br/>
      $t=T/2$
    </td>
  </tr>
</table>
图：不同时刻下浅水波算例的水深立体图
</div>

#### 5.2.2 coutour 图
<div align="center">
<table>
  <tr>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_cell_100x100_contour.png" width="400"/>
      <br/>
      $t=0$
    </td>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_T0.747541_cell_100x100_contour.png" width="400"/>
      <br/>
      $t=T/6$
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_T1.495081_cell_100x100_contour.png" width="400"/>
      <br/>
      $t=T/3$
    </td>
    <td align="center">
      <img src="2D_shallowwater_COS_DG_oscilating_lake\output\hPb_T2.242622_cell_100x100_contour.png" width="400"/>
      <br/>
      $t=T/2$
    </td>
  </tr>
</table>
图：不同时刻下浅水波算例的水深等值云图
</div>


