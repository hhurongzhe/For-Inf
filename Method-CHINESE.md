# 平面波表象下的两体矩阵元

经常遇到的是分波表象(LSJ表象)或谐振子基表象下的矩阵元，但是对于核物质的计算，经常需要平面波表象下的矩阵元作为输入:
$$
\bra{\mathbf{p}_1' s_1' m_{s1}' \mathbf{p}_2' s_2' m_{s2}'}\hat{v}\ket{\mathbf{p}_1 s_1 m_{s1} \mathbf{p}_2 s_2 m_{s2}}
$$
其中$\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_1',\mathbf{p}_2'$是两个粒子单独的初末动量；$s_1',s_2',s_1,s_2=1/2$为固定值，因此之后都隐含不写；而$m_{s1},m_{s2},m_{s1}',m_{s2}'$是自旋投影(只可能取$\pm1/2$)。

定义相对动量以及总动量:
$$
\mathbf{p}=\frac{1}{2}\left(\mathbf{p}_1-\mathbf{p}_2\right)
$$

$$
\mathbf{P}=\mathbf{p}_1+\mathbf{p}_2
$$

可以得到:
$$
\bra{\mathbf{p}_1' m_{s1}' \mathbf{p}_2' m_{s2}'}\hat{v}\ket{\mathbf{p}_1 m_{s1} \mathbf{p}_2 m_{s2}}=\bra{\mathbf{p'}\mathbf{P'} m_{s1}' m_{s2}'}\hat{v}\ket{\mathbf{p}\mathbf{P} m_{s1} m_{s2}}
$$
由于核力保持总动量以及总电荷守恒，因此:
$$
\begin{aligned}
\langle\mathbf{p'}\mathbf{P'} m_{s1}' m_{s2}'|\hat{v}|\mathbf{p}\mathbf{P} m_{s1} m_{s2}\rangle&=\delta_{T_z',T_z}\delta(\mathbf{P'}-\mathbf{P})\langle\mathbf{p'}m_{s1}' m_{s2}'|\hat{v}|\mathbf{p} m_{s1} m_{s2}\rangle\\
&=\delta_{T_z',T_z}\delta(\mathbf{P'}-\mathbf{P})\langle m_{s1}' m_{s2}'|\hat{v}(\mathbf{p'},\mathbf{p})| m_{s1} m_{s2}\rangle
\end{aligned}
$$
而上面的$\hat{v}(\mathbf{p'},\mathbf{p})$即为相互作用在动量空间、相对坐标系下的算符表达式。

例如:
$$
\hat{v}_{\mathrm{LO}}(\mathbf{p'},\mathbf{p})=\dfrac{1}{(2\pi)^3}\sqrt{\dfrac{M_N}{E_{p'}}}\Big(C_S+C_T\vec{\sigma}_1\cdot\vec{\sigma}_2-\dfrac{g_A^2}{4f_\pi^2} \mathbf{\tau}_1\cdot\mathbf{\tau}_2 \dfrac{\vec{\sigma}_1\cdot\vec{q}\;\vec{\sigma}_2\cdot\vec{q}}{q^2+m_\pi^2}\Big)\sqrt{\dfrac{M_N}{E_{p}}}\exp{\big(-(p'/\Lambda)^{2n}-(p/\Lambda)^{2n}\big)}
$$
其中$\mathbf{q}=\mathbf{p}'-\mathbf{p}$为动量转移。

因此可以直接把矩阵元的解析形式求解出来。

利用自动化分波分解(aPWD)的技巧，将相互作用进行算符分解:
$$
\hat{v}(\mathbf{p'},\mathbf{p})=\sum_{j=1}^{6}f_j\cdot w_j
$$
其中算符空间的基函数$w_j$吸收了相互作用的所有自旋依赖部分:
$$
\begin{aligned}
w_1&=1\\
w_2&=\vec{\sigma}_1 \cdot \vec{\sigma}_2\\
w_3&=i(\vec{\sigma}_1 + \vec{\sigma}_2)\cdot (\mathbf{p}\times\mathbf{p'})\\
w_4&=\vec{\sigma}_1\cdot (\mathbf{p}\times\mathbf{p'})\,\vec{\sigma}_2\cdot (\mathbf{p}\times\mathbf{p'})\\
w_5&=\vec{\sigma}_1\cdot (\mathbf{p'}+\mathbf{p})\,\vec{\sigma}_2\cdot (\mathbf{p'}+\mathbf{p})\\
w_6&=\vec{\sigma}_1\cdot (\mathbf{p'}-\mathbf{p})\,\vec{\sigma}_2\cdot (\mathbf{p'}-\mathbf{p})\\
\end{aligned}
$$
那么:
$$
\langle m_{s1}' m_{s2}'|\hat{v}(\mathbf{p'},\mathbf{p})| m_{s1} m_{s2}\rangle
=\sum_{j=1}^{6}f_j \langle m_{s1}' m_{s2}'|w_j| m_{s1} m_{s2}\rangle
$$
利用Mathematica，可以将$\langle m_{s1}' m_{s2}'|w_j| m_{s1} m_{s2}\rangle$可能的$16\times 6$种情况解析地求出来（仅仅是简单的4维矩阵运算），并且进行求和。关于$w_j$求和后只剩下16中种情况，其中一种为:
$$
\langle ++|\hat{v}(\mathbf{p'},\mathbf{p})| ++\rangle=
f_1 + f_2 + f_5 p_z'^2 + f_6 p_z'^2 + f_4 p_y'^2 p_x^2 - 2 f_4 p_x' p_y' p_x p_y + f_4 p_x'^2 p_y^2 + 2 f_5 p_z' p_z - 2 f_6 p_z' p_z + f_5 p_z^2 + f_6 p_z^2
$$
自旋投影只用正负号来进行标记。

实际上，虽然16种情况的表达式都是解析的，但是借助Mathematica的Form[]函数，可以自动生成表达式所对应的C语言代码。因此我们从来不需要真正去看每种情况的表达式是什么。我将这种方法起名为：aPWP (automated Plane-Wave Projection).

对于任何一个相互作用，只需要进行算符分解即可。

如果相互作用的算符结构十分复杂，那么进行算符分解也会是一件较为麻烦的事情 (虽然我用Mathematica也写了一个符号计算算符分解的代码)。比较幸运的是，对于手征核力而言，其算符结构是非常简单的，基本上根据表达式都可以直接写出算符的分量 (比如领头阶就只有$f_1,f_2,f_6$的成分)。

最后，数值计算时会对动量空间进行离散化，因此之前的$\delta(\mathbf{P'}-\mathbf{P})$也需要进行离散化。以均匀撒点方案为例，这会带来一个$\Delta p^3=\Big(\dfrac{2\pi}{L}\Big)^3$的因子，以及一个用来检查总动量守恒的无量纲$\delta$函数。

