#include "qpolycurve.h"

#include <QPainter>
#include <QPen>

#include "Bezier/bezier.h"

bool qPolyCurve::getDraw_control_points() const { return draw_control_points; }

void qPolyCurve::setDraw_control_points(bool value) { draw_control_points = value; }

bool qPolyCurve::getDraw_curvature_radious() const { return draw_curvature_radious; }

void qPolyCurve::setDraw_curvature_radious(bool value) { draw_curvature_radious = value; }

int qPolyCurve::type() const { return QGraphicsItem::UserType + 2; }

void qPolyCurve::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
  Q_UNUSED(option)
  Q_UNUSED(widget)

  setFlag(GraphicsItemFlag::ItemIsSelectable, true);

  painter->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform, true);
  QPen pen(Qt::darkBlue);
  pen.setStyle(isSelected() ? Qt::DashDotLine : Qt::SolidLine);
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
    painter->setPen(Qt::green);
    for (double t = 0; t <= size(); t += 1.0 / 500)
    {
      painter->setPen(QColor(static_cast<int>(std::fabs(255 * (0.5 - t / size()))),
                             static_cast<int>(255 * t / size()), static_cast<int>(255 * (1 - t / size()))));
      auto p = valueAt(t);
      auto n1 = p + normalAt(t, false) * curvatureDerivativeAt(t);
      auto n2 = p - normalAt(t, false) * curvatureDerivativeAt(t);
      painter->drawLine(QLineF(n1.x(), n1.y(), n2.x(), n2.y()));
    }
  }
}

QRectF qPolyCurve::boundingRect() const
{
  auto bbox = boundingBox(false);
  return QRectF(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(), bbox.max().y()));
}
