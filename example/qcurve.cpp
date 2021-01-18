#include "qcurve.h"

#include <QPainter>
#include <QPen>

void qCurve::setDraw_control_points(bool value) { draw_control_points = value; }

void qCurve::setDraw_curvature_radious(bool value) { draw_curvature_radious = value; }

bool qCurve::getDraw_control_points() const { return draw_control_points; }

bool qCurve::getDraw_curvature_radious() const { return draw_curvature_radious; }

bool qCurve::getLocked() const
{
    return locked;
}

void qCurve::setLocked(bool value)
{
    locked = value;
}

int qCurve::type() const { return QGraphicsItem::UserType + 1; }

void qCurve::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
  Q_UNUSED(option)
  Q_UNUSED(widget)

  setFlag(GraphicsItemFlag::ItemIsSelectable, true);

  painter->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform, true);

  QPen pen;
  pen.setStyle(isSelected() ? Qt::DashDotLine : Qt::SolidLine);
  pen.setColor(getLocked() ? Qt::red : Qt::black);
  painter->setPen(pen);
  QPainterPath curve;
  auto poly = polyline();
  curve.moveTo(poly[0].x(), poly[0].y());
  for (unsigned k = 1; k < poly.size(); k++)
    curve.lineTo(poly[k].x(), poly[k].y());
  painter->drawPath(curve);

  if (draw_control_points)
  {
    const int d = 6;
    painter->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    Bezier::PointVector points = controlPoints();
    for (unsigned k = 1; k < points.size(); k++)
    {
      painter->setPen(Qt::blue);
      painter->drawEllipse(QRectF(points[k - 1].x() - d / 2, points[k - 1].y() - d / 2, d, d));
      painter->setPen(QPen(QBrush(Qt::gray), 1, Qt::DotLine));
      painter->drawLine(QLineF(points[k - 1].x(), points[k - 1].y(), points[k].x(), points[k].y()));
    }
    painter->setPen(Qt::blue);
    painter->drawEllipse(QRectF(points.back().x() - d / 2, points.back().y() - d / 2, d, d));
  }

  if (draw_curvature_radious)
  {
    for (double t = 1.0 / 100; t <= 1.0; t += 1.0 / 200)
    {
      painter->setPen(QColor(abs(255 * (0.5 - t)), (int)(255 * t), (int)(255 * (1 - t))));
      auto p = valueAt(t);
      auto tangent = tangentAt(t);
      Bezier::Point normal(-tangent.y(), tangent.x());
      double kappa = curvatureAt(t);
      auto n1 = p + normal * kappa * 100;
      auto n2 = p - normal * kappa * 100;
      painter->drawLine(QLineF(n1.x(), n1.y(), n2.x(), n2.y()));
    }
  }
}

QRectF qCurve::boundingRect() const
{
  auto bbox = boundingBox(false);
  return QRectF(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(), bbox.max().y()));
}
