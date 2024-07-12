function [t] = chk_quantizer(r,c)
if r <= c
    return;
else 
    d = r - c;
    if c > 1
        t = (c-1) + (d>2);
    else
        t = (d>1)*(c==1);
    end
end
