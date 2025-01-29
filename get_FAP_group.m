function [flag] = get_FAP_group(FAP,mpc1,mpc2,mpc3)

flag = 4;
for stridx = 1: length(mpc1)
        if (strcmp(FAP.InteractionSummary, mpc1{stridx}))
            flag = 1;
            break
        end
end

for stridx = 1: length(mpc2)
        if (strcmp(FAP.InteractionSummary, mpc2{stridx}))
            flag = 2;
            break
        end
end


for stridx = 1: length(mpc3)
        if (strcmp(FAP.InteractionSummary, mpc3{stridx}))
            flag = 3;
            break
        end
end
return 