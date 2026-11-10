# 自由液膜破碎

二维膜面的长波模拟与热毛细波验证。当前模型推进液–气–液体系中的对称夹层厚度 $h(x,y,t)$，采用平膜迁移率、非线性范德华压力和隐式中点时间推进。这里的验证覆盖接触前的近乎平膜阶段。

## 二维无涨落正则模

在 $64\times64$ 周期网格上，取 $L_x/h_0=L_y/h_0=128$、$\Delta t^*=8$ 和吸引系数 $a=0.0084375$。初始扰动为斜向 $(m_x,m_y)=(2,2)$ 的余弦波，振幅 $10^{-3}h_0$，热噪声为零。JAM 和 Basilisk 拟合到的增长率均为 $s^*\approx1.80083\times10^{-3}$，连续长波理论给出 $1.80299\times10^{-3}$。

![二维无涨落失稳的厚度场及振幅增长](postprocess/figures/normal_mode_2d.png)

[正则模.ipynb](postprocess/正则模.ipynb) 用保存的二维厚度场和振幅历史重画上图。[tcw.ipynb](postprocess/tcw.ipynb) 比较持续加噪与首步加噪的二维壳平均谱，理论曲线采用连续波数和连续时间。首步加噪理论以实测首谱 $H_0$ 为初值。 [正则模 PDF](postprocess/figures/normal_mode_2d.pdf) · [TCW 谱 PDF](postprocess/figures/tcw_continuous.pdf)

## 仓库结构

- `jam/` 保存 Python 长波求解器。二维版本用解析雅可比求解中点方程。
- `basilisk/` 保存 Basilisk C 长波求解器和二维 FFT。二维中点版本与 JAM 接受相同的初态及噪声增量。
- `postprocess/` 保存两个可执行 notebook、紧凑的示例数据和配对计算入口。示例 TCW 数据来自 $64^2$ 网格和每种加噪方式 64 对轨迹。

## 复算

需要 NumPy、Matplotlib、Jupyter Python 内核和用于图中文字渲染的 LaTeX。使用 Basilisk `qcc` 编译二维中点程序。

```sh
python -m pip install -r requirements.txt
mkdir -p build
qcc -O3 -disable-dimensions -I"$PWD/basilisk" basilisk/gas_film_midpoint_2d.c -o build/gas_film_midpoint_2d -lm
python postprocess/run_tcw.py build/gas_film_midpoint_2d runs/tcw_512 --n 512 --replicas 30
```

$512^2$ 算例使用 $L_x/h_0=L_y/h_0=1024$，保留 $\Delta x/h_0=\Delta y/h_0=2$。在相同波数上限 $|\boldsymbol q|h_0\leq\pi/6$ 内，共有 84 个完整的非零 Fourier 壳。终点 $t^*=1280$，每条轨迹 160 步。两种协议各有 30 对轨迹，每对包含 JAM 和 Basilisk。

本机用 $512^2$ 网格实测每种协议 3 次、每次 32 步。按同一进程串行外推，持续加噪约 15 分钟，首步加噪约 8 分钟。外推覆盖推进和噪声输入准备，完整输出和谱统计另有少量开销。 [短测计时](postprocess/data/benchmark.json)
