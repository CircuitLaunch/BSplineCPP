add_library(BSpline STATIC
  BSpline.cpp BSpline.hpp
  Functor.cpp Functor.hpp
  Legendre.cpp Legendre.hpp
  Newton.cpp Newton.hpp
  Parametizer.cpp Parametizer.hpp
)

target_include_directories(BSpline PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
