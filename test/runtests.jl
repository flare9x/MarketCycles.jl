using MarketCycles
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@test 1 == 1

# Load Ehlers TS version
# TO DO
# if our output == TS output test == pass (note must round result same as TS)
