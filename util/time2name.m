function [ name ] = time2name()
%TIME2NAME Creates a string from the currrent time

time = clock;
name = [num2str(time(1)) '-' num2str(time(2)) '-' ... 
        num2str(time(3)) '-' num2str(time(4)) '-' num2str(time(5))];

end

