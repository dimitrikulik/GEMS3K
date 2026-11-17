add_executable(pstall_default_test pstall_default_test.cpp)

if(TARGET GEMS3K-static)
    target_link_libraries(pstall_default_test PRIVATE GEMS3K-static)
else()
    target_link_libraries(pstall_default_test PRIVATE GEMS3K)
endif()

target_compile_definitions(
    pstall_default_test
    PRIVATE GEMS3K_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/nodearray-gem-build/TestRW"
)

add_test(NAME pstall_default_reset COMMAND pstall_default_test)

add_executable(sorption_api_test sorption_api_test.cpp)

if(TARGET GEMS3K-static)
    target_link_libraries(sorption_api_test PRIVATE GEMS3K-static)
else()
    target_link_libraries(sorption_api_test PRIVATE GEMS3K)
endif()

target_compile_definitions(
    sorption_api_test
    PRIVATE
        GEMS3K_GIBBSITE_DATA_DIR="${CMAKE_SOURCE_DIR}/../eil_scms/tests/AluSioMSA/AluSioMSA-Gibbsite-dat"
        GEMS3K_CLAYSORP_DATA_DIR="${CMAKE_SOURCE_DIR}/../eil_scms/tests/ClaySorp/ClaySorp-Mont2tlNi01-dat"
)

add_test(NAME sorption_api_test COMMAND sorption_api_test)
