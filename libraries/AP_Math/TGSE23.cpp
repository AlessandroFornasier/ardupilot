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
*/

#pragma once

#include "TGSE23.h"

const Vector9F Vector9F::operator+(const Vector9F &other) const
{
  Vector9F result;
  result.w = w + other.w;
  result.x = x + other.x;
  result.y = y + other.y;
  return result;
}

const Vector9F Vector9F::operator+=(const Vector9F &other)
{
  return *this = *this + other;
}

const Vector9F Vector9F::operator-(const Vector9F &other) const
{
  Vector9F result;
  result.w = w - other.w;
  result.x = x - other.x;
  result.y = y - other.y;
  return result;
}

const Vector9F Vector9F::operator-=(const Vector9F &other)
{
  return *this = *this - other;
}

const TriangularMatrix9F TriangularMatrix9F::operator*(const TriangularMatrix9F &other) const
{
  TriangularMatrix9F result;
  result.A1 = A1 * other.A1;
  result.A2 = A2 * other.A1 + A1 * other.A2;
  result.A3 = A3 * other.A1 + A1 * other.A3;
  return result;
}

const TriangularMatrix9F TriangularMatrix9F::operator*=(const TriangularMatrix9F &other)
{
  return *this = *this * other;
}

const TriangularMatrix9F TriangularMatrix9F::operator*(const ftype &other) const
{
  TriangularMatrix9F result;
  result.A1 = A1 * other;
  result.A2 = A2 * other;
  result.A3 = A3 * other;
  return result;
}

const TriangularMatrix9F TriangularMatrix9F::operator*=(const ftype &other)
{
  return *this = *this * other;
}

const TriangularMatrix9F TriangularMatrix9F::operator+(const TriangularMatrix9F &other) const
{
  TriangularMatrix9F result;
  result.A1 = A1 + other.A1;
  result.A2 = A2 + other.A2;
  result.A3 = A3 + other.A3;
  return result;
}

const TriangularMatrix9F TriangularMatrix9F::operator+=(const TriangularMatrix9F &other)
{
  return *this = *this + other;
}

const TriangularMatrix9F TriangularMatrix9F::operator-(const TriangularMatrix9F &other) const
{
  TriangularMatrix9F result;
  result.A1 = A1 - other.A1;
  result.A2 = A2 - other.A2;
  result.A3 = A3 - other.A3;
  return result;
}

const TriangularMatrix9F TriangularMatrix9F::operator-=(const TriangularMatrix9F &other)
{
  return *this = *this - other;
}

const Vector9F TriangularMatrix9F::operator*(const Vector9F &other) const
{
  Vector9F result;
  result.w = A1 * other.w;
  result.x = A2 * other.w + A1 * other.x;
  result.y = A3 * other.w + A1 * other.y;
  return result;
}

const SE23 SE23::identity()
{
  Matrix3F R;
  R.identity();
  return SE23(R, Vector3F(), Vector3F());
}

const SE23 SE23::exponential(const Vector9F &u)
{
  const ftype ang = u.w.length();

  const Matrix3F I = Matrix3F().identity();
  const Matrix3F W = Matrix3F::skew_symmetric(u.w);

  const ftype a, b, c;
  if (ang < _eps)
  {
    a = 1.0f;
    b = 0.5f;
    c = 1.0f / 6.0f;
  }
  else
  {
    a = sinF(ang) / ang;
    b = (1.0f - cosF(ang)) / (ang * ang);
    c = (1.0f - a) / (ang * ang);
  }

  const Matrix3F R = I + a * W + b * W * W;
  const Matrix3F JL = I + b * W + c * W * W;
  const Vector3F W1 = JL * u.x;
  const Vector3F W2 = JL * u.y;

  return SE23(R, W1, W2);
}

const SE23 SE23::operator*(const SE23 &other) const
{
  const Matrix3F R = _R * other._R;
  const Vector3F W1 = _R * other._W1 + _W1;
  const Vector3F W2 = _R * other._W2 + _W2;
  return SE23(R, W1, W2);
}

const SE23 SE23::inverse() const
{
  const Matrix3F RT = _R.transpose();
  const Vector3F W1 = -RT * _W1;
  const Vector3F W2 = -RT * _W2;
  return SE23(R, W1, W2);
}

const TriangularMatrix9F SE23::JL(const Vector9F &u)
{
  TriangularMatrix9F J;
  J.A1.identity();

  const ftype ang = u.w.length();
  const Matrix9F W = SE23::adjointMatrix(u);

  if (ang < _eps)
  {
    return J + 0.5f * W + (1.0f / 6.0f) * W * W;
  }

  const ftype a = sinF(ang) / ang;
  const ftype b = (1.0f - cosF(ang)) / (ang * ang);
  const ftype c = (1.0f - a) / (ang * ang);

  J.A1 += b * W + c * W * W;
  J.A2 = SE23::Q(u.w, u.x);
  J.A3 = SE23::Q(u.w, u.y);

  return J;
}

const TriangularMatrix9F SE23::Adjoint() const
{
  TriangularMatrix9F Ad;
  Ad.A1 = _R;
  Ad.A2 = Matrix3F::skew_symmetric(_W1) * _R;
  Ad.A3 = Matrix3F::skew_symmetric(_W2) * _R;
  return Ad;
}

const TriangularMatrix9F SE23::inverseAdjoint() const
{
  TriangularMatrix9F Adi;
  Adi.A1 = _R.transpose();
  Adi.A2 = -Ad.A1 * Matrix3F::skew_symmetric(_W1);
  Adi.A3 = -Ad.A1 * Matrix3F::skew_symmetric(_W2);
  return Adi;
}

const TriangularMatrix9F SE23::adjointMatrix(const Vector9F &u)
{
  TriangularMatrix9F adj;
  adj.A1 = Matrix3F::skew_symmetric(u.w);
  adj.A2 = Matrix3F::skew_symmetric(u.x);
  adj.A3 = Matrix3F::skew_symmetric(u.y);
  return adj;
}

const Matrix3F SE23::Q(const Vector3F &u1, const Vector3F &u2)
{
  Matrix3F p = Matrix3F::skew_symmetric(w);
  Matrix3F r = Matrix3F::skew_symmetric(v);

  const ftype ang = u.w.length();
  const ftype s = sinF(ang);
  const ftype c = cosF(ang);

  const ftype ang_p2 = powF(ang, 2);
  const ftype ang_p3 = ang_p2 * ang;
  const ftype ang_p4 = ang_p3 * ang;
  const ftype ang_p5 = ang_p4 * ang;

  const ftype c1 = (ang - s) / ang_p3;
  const ftype c2 = (0.5 * ang_p2 + c - 1.0f) / ang_p4;
  const ftype c3 = (ang * (1.0f + 0.5f * c) - 1.5f * s) / ang_p5;

  Matrix3F m1 = p * r + r * p + p * r * p;
  Matrix3F m2 = p * p * r + r * p * p - 3.0 * p * r * p;
  Matrix3F m3 = p * r * p * p + p * p * r * p;

  return 0.5 * r + c1 * m1 + c2 * m2 + c3 * m3;
}

const TGSE23 TGSE23::identity()
{
  return TGSE23(SE23::identity(), Vector9F());
}

const TGSE23 TGSE23::exponential(const Vector9F &u1, const Vector9F &u2)
{
  return TGSE23(SE23::exponential(u1), SE23::JL(u1) * u2);
}

const TGSE23 TGSE23::operator*(const TGSE23 &other) const
{
  return TGSE23(_D * other._D, _d + _D.Adjoint() * other._d);
}

const TGSE23 TGSE23::inverse() const
{
  return TGSE23(_D.inv(), -_D.inverseAdjoint() * _d);
}