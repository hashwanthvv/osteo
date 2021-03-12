%% Modeling of Osteosarcoma and mutations using differential equations
% Author: Haswanth Vundavilli
% Email: hashwanthvv@tamu.edu

  
%% Differential equations

function auc = osteosarcoma(mutation, x1, x2, x3, x4, x5, x6, x7)

   syms g1(t) g2(t) g3(t) g4(t) g5(t) g6(t) g7(t) g8(t) g9(t) g10(t) g11(t) ...
   g12(t) g13(t) g14(t) g15(t) g16(t) g17(t) g18(t) g19(t) g20(t) g21(t) ...
   g22(t) g23(t) g24(t) g25(t) g26(t) g27(t) g28(t) g29(t) g30(t) output(t)
   
    % Drugs
    Cryptotanshinone = x1; Temsirolimus = x2; LY294002 = x3; Lapatinib = x4; ...
    HO3867 = x5; NT157 = x6; U0126 = x7;
    
    k = 0.05;
    
    if (mutation == 1)
        ode1 = diff(g1) == (k)*g1; 
        cond1 = g1(0) == 1;
    else
        ode1 = diff(g1) == 0; 
        cond1 = g1(0) == 0;
    end
    
    if (mutation == 2)
        ode2 = diff(g2) == (k)*g2; 
        cond2 = g2(0) == 1;
    else
        ode2 = diff(g2) == (k)*(g1-g2); 
        cond2 = g2(0) == 0;
    end
    
    if (mutation == 3)
        ode3 = diff(g3) == (k)*g3; 
        cond3 = g3(0) == 1;
    else
        ode3 = diff(g3) == (k)*(g2-2*g3); 
        cond3 = g3(0) == 0;
    end
    
    if (mutation == 4)
        ode4 = diff(g4) == (k)*(g4)*(1-Cryptotanshinone)*(1-HO3867); 
        cond4 = g4(0) == (1-Cryptotanshinone)*(1-HO3867);
    else
        ode4 = diff(g4) == (k)*(g3-5*g4)*(1-Cryptotanshinone)*(1-HO3867); 
        cond4 = g4(0) == 0;
    end
    
    if (mutation == 5)
        ode5 = diff(g5) == (k)*g5; 
        cond5 = g5(0) == 1;
    else
        ode5 = diff(g5) == 0; 
        cond5 = g5(0) == 0;
    end
    
    if (mutation == 6)
        ode6 = diff(g6) == (k)*g6; 
        cond6 = g6(0) == 1;
    else
        ode6 = diff(g6) == (k)*(g5-g6); 
        cond6 = g6(0) == 0;
    end
    
    if (mutation == 7)
        ode7 = diff(g7) == (k)*g7*(1-NT157); 
        cond7 = g7(0) == 1-NT157;
    else
        ode7 = diff(g7) == (k)*(g6-g7)*(1-NT157); 
        cond7 = g7(0) == 0;
    end
    
    if (mutation == 8)
        ode8 = diff(g8) == (k)*(g8)*(1-LY294002); 
        cond8 = g8(0) == 1-LY294002;
    else
        ode8 = diff(g8) == (k)*(g3+g7+g12+g17-g8)*(1-LY294002); 
        cond8 = g8(0) == 0;
    end
        
    if (mutation == 9)
        ode9 = diff(g9) == (k)*(g9); 
        cond9 = g9(0) == 1;
    else
        ode9 = diff(g9) == (k)*(g8-g9); 
        cond9 = g9(0) == 0;
    end
    
    if (mutation == 10)
        ode10 = diff(g10) == (k)*(g10); 
        cond10 = g10(0) == 1;
    else
        ode10 = diff(g10) == (k)*(g9-3*g10); 
        cond10 = g10(0) == 0;
    end
    
    if (mutation == 11)
        ode11 = diff(g11) == (k)*(g11); 
        cond11 = g11(0) == 1;
    else
        ode11 = diff(g11) == 0; 
        cond11 = g11(0) == 0;
    end
    
    if (mutation == 12)
        ode12 = diff(g12) == (k)*(g12)*(1-Lapatinib); 
        cond12 = g12(0) == 1-Lapatinib;
    else
        ode12 = diff(g12) == (k)*(g11-2*g12)*(1-Lapatinib); 
        cond12 = g12(0) == 0;
    end
    
    if (mutation == 13)
        ode13 = diff(g13) == (k)*(g13); 
        cond13 = g13(0) == 1;
    else
        ode13 = diff(g13) == 0; 
        cond13 = g13(0) == 0;
    end
    
    if (mutation == 14)
        ode14 = diff(g14) == (k)*(g14); 
        cond14 = g14(0) == 1;
    else
        ode14 = diff(g14) == (k)*(g13-g14); 
        cond14 = g14(0) == 0;
    end
    
    if (mutation == 15)
        ode15 = diff(g15) == (k)*(g15); 
        cond15 = g15(0) == 1;
    else
        ode15 = diff(g15) == (k)*(g12+g14-g15); 
        cond15 = g15(0) == 0;
    end
        
    if (mutation == 16)
        ode16 = diff(g16) == (k)*(g16); 
        cond16 = g16(0) == 1;
    else
        ode16 = diff(g16) == (k)*(g15-g16); 
        cond16 = g16(0) == 0;
    end
    
    if (mutation == 17)
        ode17 = diff(g17) == (k)*(g17); 
        cond17 = g17(0) == 1;
    else
        ode17 = diff(g17) == (k)*(g16-2*g17); 
        cond17 = g17(0) == 0;
    end
    
    if (mutation == 18)
        ode18 = diff(g18) == (k)*(g18); 
        cond18 = g18(0) == 1;
    else
        ode18 = diff(g18) == (k)*(g17-g18); 
        cond18 = g18(0) == 0;
    end
    
    if (mutation == 19)
        ode19 = diff(g19) == (k)*(g19)*(1-U0126); 
        cond19 = g19(0) == 1-U0126;
    else
        ode19 = diff(g19) == (k)*(g18-g19)*(1-U0126); 
        cond19 = g19(0) == 0;
    end
    
    if (mutation == 20)
        ode20 = diff(g20) == (k)*(g20)*(1-Cryptotanshinone); 
        cond20 = g20(0) == 1-Cryptotanshinone;
    else
        ode20 = diff(g20) == (k)*(g19+g10-3*g20)*(1-Cryptotanshinone); 
        cond20 = g20(0) == 0;
    end
    
    if (mutation == 21)
        ode21 = diff(g21) == (k)*(g21)*(1-Temsirolimus); 
        cond21 = g21(0) == 1-Temsirolimus;
    else
        ode21 = diff(g21) == (k)*(g10-g21)*(1-Temsirolimus); 
        cond21 = g21(0) == 0;
    end
    
    if (mutation == 22)
        ode22 = diff(g22) == (k)*(g22); 
        cond22 = g22(0) == 1;
    else
        ode22 = diff(g22) == (k)*(g10-g22); 
        cond22 = g22(0) == 0;
    end
        
    if (mutation == 23)
        ode23 = diff(g23) == (k)*(g23); 
        cond23 = g23(0) == 1;
    else
        ode23 = diff(g23) == (k)*(g20+g21-g23); 
        cond23 = g23(0) == 0;
    end
    
     
   odes = [ode1; ode2; ode3; ode4; ode5; ode6; ode7; ode8; ode9; ode10; ...
        ode11; ode12; ode13; ode14; ode15; ode16; ode17; ode18; ode19; ...
        ode20; ode21; ode22; ode23];
    
   conds = [cond1; cond2; cond3; cond4; cond5; cond6; cond7; cond8; cond9; ...
       cond10; cond11; cond12; cond13; cond14; cond15; cond16; cond17;  ...
       cond18; cond19; cond20; cond21; cond22; cond23];   
   
    S = dsolve(odes,conds);

    g_24(t) = S.g4;
    g_25(t) = (S.g4)*(S.g10);
    g_26(t) = (S.g4)*(S.g10);
    g_27(t) = (S.g4)*(S.g22);
    g_28(t) = S.g23;
    g_29(t) = (S.g4)+(S.g20);
    g_30(t) = S.g20;
   
    
    output(t) = g_24(t) + g_25(t) + g_26(t) + g_27(t) + g_28(t) + g_29(t) + g_30(t);
    area = int(output(t),[0 100]);
    auc = double(area);
    
end
