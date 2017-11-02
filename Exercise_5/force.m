function F = force(x,h)
    F = -6.23347*exp(-3.125*(-5 + x)^2)*(-5 + x)*h;