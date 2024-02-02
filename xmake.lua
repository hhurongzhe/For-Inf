add_rules("mode.release")

target("apwp")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("aggressive")
    add_files("src/aPWP.cpp")