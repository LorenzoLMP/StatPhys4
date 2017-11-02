function V = potential(x, h)
    V = h*normpdf(x,5,0.4);
    