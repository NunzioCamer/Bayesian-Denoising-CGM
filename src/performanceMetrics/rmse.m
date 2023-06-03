function [RMSE] = rmse(y1,y2)
%y1 e y2 input signals to compare

RMSE = sqrt(mean((y1-y2).^2));

end
