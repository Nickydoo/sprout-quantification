function [xout, yout] = untangle(xin,yin)

xout = xin(1);
yout = yin(1);

xin(1) = [];
yin(1) = [];
while not(length(xin)==0)            
    d = (xout(end)-xin).^2 + (yout(end)-yin).^2;
    [minD,iMin] = min(d);
    if minD>2
        xout = [xout;NaN];
        yout = [yout;NaN];
        xout = [xout;xin(iMin)];
        yout = [yout;yin(iMin)];
    else
        xout = [xout;xin(iMin)];
        yout = [yout;yin(iMin)];
    end
    xin(iMin) = [];
    yin(iMin) = [];
end
