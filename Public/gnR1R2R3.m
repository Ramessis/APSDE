function [r1,r2,r3] = gnR1R2R3(NP1, r0)
    if NP1 > 4
        NP0=length(r0);
        for i = 1:NP0
            temp = randperm(NP1);
            r1r2r3=temp(1:3);
            r1r2r3(find(r1r2r3 == r0(i))) = temp(4);
            r1(i)=r1r2r3(1);
            r2(i)=r1r2r3(2);
            r3(i)=r1r2r3(3);
        end
    end
            
end 