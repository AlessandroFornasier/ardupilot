/*
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2024 Alessandro Fornasier (alessandro.fornasier@ieee.org),
    all rights reserved.

    Implementation of the Extended Special Euclidean group SE2(3).
*/

#pragma once

#include "AP_Math.h"
#include "ftype.h"
#include "vector3.h"
#include "vectorN.h"
#include "matrix3.h"
#include "matrixN.h"

struct Vector9F
{
  const Vector9F operator+(const Vector9F &other) const;
  const Vector9F operator+=(const Vector9F &other);
  const Vector9F operator-(const Vector9F &other) const;
  const Vector9F operator-=(const Vector9F &other);

  Vector3F w;
  Vector3F x;
  Vector3F y;
};

// Block triangular matrix for SE23 Adjoint, SE23 Left Jacobian and se23 adjoint
struct TriangularMatrix9F
{
  // this *+- other
  const TriangularMatrix9F operator*(const TriangularMatrix9F &other) const;
  const TriangularMatrix9F operator*=(const TriangularMatrix9F &other);
  const TriangularMatrix9F operator*(const ftype &other) const;
  const TriangularMatrix9F operator*=(const ftype &other);
  const TriangularMatrix9F operator+(const TriangularMatrix9F &other) const;
  const TriangularMatrix9F operator+=(const TriangularMatrix9F &other);
  const TriangularMatrix9F operator-(const TriangularMatrix9F &other) const;
  const TriangularMatrix9F operator-=(const TriangularMatrix9F &other);

  // this * vector
  const Vector9F operator*(const Vector9F &other) const;

  Matrix3F A1;
  Matrix3F A2;
  Matrix3F A3;
};

class SE23
{
public:
  // default constructor
  SE23() = default;

  // Element contructor
  SE23(const Matrix3F &R, const Vector3F &W1, const Vector3F &W2) : _R(R), _W1(W1), _W2(W2){};

  // Identity element
  static const SE23 identity();

  // Exponential map for SE23
  static const SE23 exponential(const Vector9F &u);

  // Group operation this * other
  const SE23 operator*(const SE23 &other) const;

  // Inverse
  const SE23 inverse() const;

  // SE23 Adjoint matrices
  const TriangularMatrix9F Adjoint() const;
  const TriangularMatrix9F inverseAdjoint() const;

  // se23 adjoint matrix
  static const TriangularMatrix9F adjointMatrix(const Vector9F &u);

  // SE23 Left Jacobian
  static const TriangularMatrix9F JL(const Vector9F &u);

private:
  // Q matrix for SE23 Left jacobian
  static const Matrix3F Q(const Vector3F &u1, const Vector3F &u2);

  Matrix3F _R;
  Vector3F _W1;
  Vector3F _W2;

  static constexpr ftype _eps = std::is_same_v<ftype, float> ? 1.0e-6f : 1.0e-9f;
};

class TGSE23
{

public:
  // default constructor
  TGSE23() = default;

  // Element contructor
  TGSE23(const SE23 &D, const Vector9F &d) : _D(D), _d(d){};

  // Identity element
  static const TGSE23 identity();

  // Exponential map for SE23
  static const TGSE23 exponential(const Vector9F &u1, const Vector9F &u2);

  // Group operation
  const TGSE23 operator*(const TGSE23 &other) const;

  // Inverse
  const TGSE23 inverse() const;

private:
  SE23 _D;
  Vector9F _d;

  static constexpr ftype _eps = std::is_same_v<ftype, float> ? 1.0e-6f : 1.0e-9f;
};