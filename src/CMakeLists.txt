add_library(BSpline STATIC
  BSpline.cpp BSpline.hpp
  Legendre.cpp Legendre.hpp
  Newton.cpp Newton.hpp
  Parametizer.cpp Parametizer.hpp Functor.hpp
)

target_include_directories(BSpline PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
