---
title: The Magnus Effect
math: true
date : 2024-01-24
draft : true
---

# Magnus effect 

This article aims to further strengthen the understanding of the Magnus
effect by finding a direct correlation between the Bernoulli Fluid
Equation and the Magnus force. The Magnus force is force that affects
\"rapidly\" spinning bodies in the direction perpendicular to their
movement. Now here comes the construction from the Bernoulli equation.
```katex
\int_{0}^{m}dx
```

## Analytical approach

Suppose a cylinder moves forward so the \"faster\" moving surface
\"pulls\" the air up and in turn gets pulled down due to newton's Third
law.

## Scenario

There is a cylinder moving in space in the direction $$+\hat{x}$$. Now the
cylinder is given an angular velocity $$\omega$$ which points in the
$$+\hat{k}$$ direction.Therefore the ball spins counter clockwise.

## Derivation

We start of with the Bernoulli equation from fluid dynamics. In this
case the terms are as follows: $$P_0$$ refers to the ambient pressure,
$$\frac{1}{2}\rho v^2$$ refers to the dynamic pressure, the last term can
be ignored here as the cylinder used for measuring the Magnus effect can
be taken to be arbitrarily small causing the height difference to be
negligible. $$P_0 + \frac{1}{2}\rho v^2 +\rho_0 gh = \kappa $$ The
ambient pressure change essentially gives rise to the change in the
pressure between the high and low spinning sides of the cylinder causing
the Magnus effect. So here the term $v$ represents the fluid velocity.
here the fluid velocity is the velocity is the velocity of a small
elemental part of the fluid. Now let us define two pressures of the
cylinder in question. Let the top pressure be denoted by $P_t$ and the
bottom surface pressure be denoted by $P_b$. Now let the ball move from
left to right in this question (suppose $$\hat{x}$$ with velocity u).
$$
    P_t = \kappa - \frac{1}{2}\rho(v-r\omega)^2 $$ $$
    P_b = \kappa - \frac{1}{2}\rho(v+r\omega)^2
$ To understand the above equations let us enter the COM
frame of the object(only translation frame not the spinning frame). This
essentially allows us to set the velocity of air in direction
$-\hat{x}$. So now we need the fluid velocity at the top of the cylinder
and the bottom of the cylinder. The fluid velocity is the relative
velocity of the air with respect to the top surface of the cylinder. At
the top part of the cylinder we have $v_r = v-r\omega$ as the ball is
spinning in direction opposite to the velocity at that point. Similarly
the expression for the bottom velocity can also be defined.
$$\begin{aligned}
    \Delta P &= |-\frac{1}{2}\rho(v-r\omega)^2 +\frac{1}{2}\rho(v+r\omega)^2|\\
    \Delta P &= |2\rho v r \omega|
\end{aligned}$$ This $\Delta P$ is in the net direction downward i.e. in
the direction $\omega \times v$. This net pressure change gives us an
important fact that there will be a net force downward causing the ball
to 'curl'.

## Conclusion

So essentially we can define the Magnus force as: $$\begin{aligned}
\fbox{$\vb{F_{magnus}} = S(v)\vb*{\omega\times v}$}
\end{aligned}$$ Now there is a lack of published formula on the Magnus
Effect because there are multiple formulae that yield different results.
We can assume irrotational flow which gives us the above definition of
the Magnus effect. But , without assuming the irrotationality of the
flow we can assume that the liquid is being pulled away due to the
centrifugal force and thus needs an added pressure gradient to balance
out the centrifugal force experienced by the liquid. Further on that
later.

## Post Notes

The defined formula for the Magnus effect has not been published yet as
the defined formulae come out to be different for different models of
the fluid flow taken and the assumptions made in the liquid's movement.
The general formula written above is True but the formula for $S(v)$ is
not defined properly.

## Extra Reading

1.  Irrotational Flow in Fluids(Richard Fitzpatrick)

2.  Potential Flow around a cylinder

3.  David Tong (Fluid Dynamics)

4.  Flow past a circular cylinder
