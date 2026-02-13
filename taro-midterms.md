Standard basis: $\{1,t,t^2,\dots,t^d\};$ Shifted Basis: $S_c=\{1,t-c,(t-c)^2,\dots,(t-c)^d\}$
Vandermonde Basis: $V(t_0,t_1,\dots,t_d)=\{(t-t_0)^d,(t-t_1)^d,(t-t_2)^d,\dots,(t-t_d)^d\}$
Bernstein Basis: $\{B_0^d(t),B_1^d(t),\dots,B_d^d(t)\};B_i^d(t)=\binom di(1-t)^{d-i}\cdot t^i;\binom di=\frac{d!}{i!(d-i)!}$
Cumulative Bernstein: $\{C_0^d(t),C_1^d(t),\dots,C_d^d(t)\};C_{i}^{d}(t)=\sum_{j=i}^{d}B_{j}^{d}(t)$
BB/CB-form: $\gamma(t)=\sum^d_{i=0}B^d_i(t)P_i=P_0+\sum^d_{i=1}(P_i-P_{i-1})C^d_i(t)$
NLI form: $\gamma_d(t)=(1-t)_{\gamma[P_0,P_1,\dots,P_{d-1}]}(t)+t_{\gamma[P_1,P_2,\dots,P_d]}(t)$
Midpoint: get two edges of NLI triangle, recursively apply to them
$C^d_0=\sum^d_{j=0}B^d_j(t)=1;B^d_i(t)=B^d_{d-i}(1-t);B^d_i(t)>0,\forall t\in(0,1)$
$\binom nk=\binom{d-1}{i-1}+\binom{d-1}{i};\frac{d}{dt}C^d_i(t)=dB^{d−1}_{i−1}(t)$
$B^d_i(t)=tB^{d-1}_{i-1}(t)+(1-t)B^{d-1}_i(t);\frac{d}{dt}B_{i}^{d}(t)=d\left(B_{i-1}^{d-1}(t)-B_{i}^{d-1}(t)\right)$
Confluent Vandermonde Determinants:
$D(abc)=\left|\begin{matrix}1&a&a^2\\1&b&b^2\\1&c&c^2\end{matrix}\right|=(b-a)(c-a)(c-b);D(a^2b)=\left|\begin{matrix}1&a&a^2\\0&1&2a\\1&b&b^2\end{matrix}\right|=(b-a)^2$
In general: $D(a_1^{e_1}a_2^{e_2}\dots a_n^{e_n})=\prod_{1\le i\le j \le n}(a_i-a_j)^{e_i\cdot e_j}\cdot \prod^n_{k=1}(e_k-1)!!$
and $k!!=\prod_{i=0}^k(i!);$ E.g., $2!!=2,3!!=12,4!!=24\cdot12=288$
Lagrange Poly: $P_n(t)=\sum^n_{k=0}l^n_k(t)\cdot f(t_k);l^n_k(t)=\prod^n_{i\ne k}\frac{t-t_i}{t_k-t_i}$
Newton Basis: $\{N_0(t),N_1(t),\dots,N_n(t)\};N_0(t)=1;N_i(t)=N_{i-1}(t)\cdot(t-t_{i-1})$
Newton Form: $p(t)=\sum_{i=0}^d[t_0,\dots,t_i]_gN_i(t)$
$[t_i]_g=g(t_i);[t_0,t_1,\dots,t_k]_g=\frac{[t_1,\dots,t_k]_g-[t_0,\dots,t_{k-1}]_g}{t_k-t_0}$
$[t_i,t_{i+1},\dots,t_{i+k}]_g=\frac{[t_{i+1},\dots,t_{i+k}]_g-[t_{i},\dots,t_{i+k-1}]_g}{t_{i+k}-t_i}$
$[t_i,t_{i+1},\dots,t_{i+k}]_g=g^{(k)}(t)$ (for repeated roots with multiplicity $\le k$)
Ordered $k$-tuples ($P^k_d$): cartesian product of ${P_d \times P_d \times \dots \times P_d} \text{ k times}$
$P_d^k$ standard basis: the union set of standard bases in each dimension
Shifted and Truncated Power: $(x-c)^k;(x-c)^k_+=\begin{cases}(x-c)^k,&x>c\\0,&x\le c\end{cases}$
${P_{d,r}^k}_{[0,u_1,.,u_k]}=\{1,t,.,t^d,(t-u_1)^{r+1}_+,.,(t-u_1)^d_+,.,(t-u_{k-1})^{r+1}_+,.,(t-u_{k-1})^d_+\}$
Order of Continuity $r$: $f(t)\in C^{r}[a,b]:f(t),f'(t),\dots,f^{(r)}(t)$ are all continuous in $[a,b]$
Piecewise function in Truncated Power form: get each interval $[c_i,c_{i+1})$ and poly $P_i(t);$
$1_{[c,\infty)}(t)=(x-c)^0_1;1_{(-\infty,c)}(t)=1-(x-c)^0_+;P(t)=\sum^n_{i=0}a_ix^i$
$1_{[c_i,c_{i+1})}(t)=(t-c_i)^0_+-(t-c_{i+1})^0_+;f(t)=\sum^n_{i=0}1_{[c_i,c_{i+1})}(t)\cdot P_i(t)$
$x^n=\sum^n_{j=0}\binom njc^{n-j}(x-c)^j;x^n(x-c)^k_+=\sum^n_{j=0}\binom njc^{n-j}(x-c)^{j+k}_+$
Correspondence between $P^k_d$ and ${P_d^k}_{[u_0,\dots,u_k]}$: map each member in $k$-tuple with an interval
$[t_i,t_{i+1},\dots,t_{i+k}]_f=\sum_{r=i}^{i+k}([t_i,\dots,t_r]_g)([t_r,\dots,t_{i+k}]_h)$ for $f(t)=g(t)\cdot h(t)$
If no solution (error, not enough degree), approx: $V\vec c\approx\vec y;V^TV\vec c=V^T\vec y$
