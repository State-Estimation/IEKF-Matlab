% Function to shift the angle
function shifted = shift_to_final(goal, current)
    shifted=current;
    len=length(goal);
    for i=1:len
        if current(i) < (goal(i) - pi)
            current(i) = current(i) + 2*pi;
        elseif current(i) > (goal(i) + pi)
            current(i) = current(i) - 2*pi;
        end
        shifted(i) = current(i);
    end
end