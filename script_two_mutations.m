n = 2; % The number of drug combinations

% Overall measure for single mutation network
A = zeros(23,23,128);

% "i" and "j" iterates over faults, and d1, d2, d3, d4, d5 iterate over drugs
for i = 1:23
    for j = 1:23
        for d1 = 0:1
            for d2 = 0:1
                for d3 = 0:1
                    for d4 = 0:1
                        for d5 = 0:1
                            for d6 = 0:1
                                for d7 = 0:1
                                    m = 64*d1 + 32*d2 + 16*d3 + 8*d4 + 4*d5 + 2*d6 + d7 + 1;
                                    if (d1+d2+d3+d4+d5+d6+d7 <= n) 
                                       A(i,j,m) = osteosarcoma_two_mutations(i,j,d1,d2,d3,d4,d5,d6,d7);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
