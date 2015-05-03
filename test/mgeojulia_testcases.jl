################################################################################
#                          Test cases for MGEOjulia.
################################################################################
#
# This file defines the test cases for MGEOjulia.
# The scenarios were obtained from:
#
#     Galski, R. L (2006). Desenvolvimento de Versões Aprimoradas, Híbridas,
#     Paralela e Multiobjetivo do Método da Otimização Extrema Generalizada e
#     Sua Aplicação no Projeto de Sistemas Espaciais. Ph.D. Thesis. Instituto
#     Nacional de Pesquisas Espaciais: Brazil.
#
################################################################################

using MGEOjulia

################################################################################
# Objective functions
################################################################################

function f_obj1(vars::Array{Float64,1})
    x1 = vars[1]

    f = [exp(-x1) + 1.4*exp(-x1*x1);
         exp(+x1) + 1.4*exp(-x1*x1)]

    (true, f)
end

function f_obj2(vars::Array{Float64,1})
    x1 = vars[1]

    f = [x1*x1;
         (x1-2.0)*(x1-2.0)]

    (true, f)
end

function f_obj3(vars::Array{Float64,1})
    x1 = vars[1]

    f = Array(Float64, 2)

    if (x1 <= 1.0)
        f[1] = -x1
    elseif ( (1.0 < x1) && (x1 <= 3.0) )
        f[1] = -2.0+x1
    elseif ( (3.0 < x1) && (x1 <= 4.0) )
        f[1] = 4.0-x1
    else
        f[1] = -4.0+x1
    end
        
    f[2] = (x1-5)*(x1-5);

    (true, f)
end

function f_obj4(vars::Array{Float64,1})
    x = vars[1]
    y = vars[2]

    f = [x*x + y*y;
         (x-1.5)*(x-1.5)+(y-3.5)*(y-3.5)]

    (true, f)
end

################################################################################
# Test cases
################################################################################

testcase = 2

if (testcase == 1)

    designVars = confDesignVars(1, 16, -10.0, +10.0)
    mgeoData = MGEOStructure(2, 0.5, 8000, 50, designVars, 1e-10)
    @time pf = MGEOrun(mgeoData, f_obj1, true)
    sortParetoFrontier(pf, 1)
    savePFtoCSV("testcase_01.csv",designVars, pf)
    
elseif (testcase == 2)

    designVars = confDesignVars(1, 16, -10.0, +10.0)
    mgeoData = MGEOStructure(2, 0.5, 8000, 50, designVars, 1e-10)
    @time pf = MGEOrun(mgeoData, f_obj2, true)
    sortParetoFrontier(pf, 1)
    savePFtoCSV("testcase_02.csv",designVars, pf)
    
elseif (testcase == 3)

    designVars = confDesignVars(1, 16, -10.0, +10.0)
    mgeoData = MGEOStructure(2, 0.5, 8000, 50, designVars, 1e-10)
    @time pf = MGEOrun(mgeoData, f_obj3, true)
    sortParetoFrontier(pf, 1)
    savePFtoCSV("testcase_03.csv",designVars, pf)
    
elseif (testcase == 4)

    designVars = confDesignVars(2, 8, -10.0, +10.0)
    mgeoData = MGEOStructure(2, 0.5, 8000, 50, designVars, 1e-10)
    @time pf = MGEOrun(mgeoData, f_obj4, true)
    sortParetoFrontier(pf, 1)
    savePFtoCSV("testcase_04.csv",designVars, pf)
    
else
    println("Invalid test case.")
end 
