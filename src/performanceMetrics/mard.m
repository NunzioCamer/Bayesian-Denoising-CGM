function [MARD] = mard(y1,y2)
%y1: filtered signal
%y2: ground-truth signal

MARD = mean(abs(y1-y2)./y2*100);

end

