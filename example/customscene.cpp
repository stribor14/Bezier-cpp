#include "customscene.h"

#include <QGraphicsSceneMouseEvent>
#include <QKeyEvent>
#include <QMessageBox>

#define is_curve (curve->type() == QGraphicsItem::UserType + 1)
#define is_curve_item(x) (x->type() == QGraphicsItem::UserType + 1)
#define is_poly (curve->type() == QGraphicsItem::UserType + 2)

#define c_curve (static_cast<qCurve*>(curve))
#define c_curve_item(x) (static_cast<qCurve*>(x))
#define c_poly (static_cast<qPolyCurve*>(curve))

void CustomScene::drawForeground(QPainter* painter, const QRectF& rect)
{
  Q_UNUSED(rect)

  if (draw_box_)
  {
    painter->setPen(Qt::blue);
    for (auto&& curve : items())
    {
      Bezier::BoundingBox bbox;
      if (is_curve)
        bbox = c_curve->boundingBox(true);
      if (is_poly)
        bbox = c_poly->boundingBox(true);
      painter->drawRect(bbox.min().x(), bbox.min().y(), bbox.max().x() - bbox.min().x(),
                        bbox.max().y() - bbox.min().y());
    }
  }

  if (draw_inter_)
  {
    painter->setPen(Qt::red);
    painter->setBrush(QBrush(Qt::red, Qt::SolidPattern));
    for (int k = 0; k < items().size(); k++)
      for (int i = k; i < items().size(); i++)
      {
        Bezier::PointVector inter;
        if (items()[k]->type() == QGraphicsItem::UserType + 1 && items()[i]->type() == QGraphicsItem::UserType + 1)
          inter = static_cast<qCurve*>(items()[i])->pointsOfIntersection(*static_cast<qCurve*>(items()[k]));
        if (items()[k]->type() == QGraphicsItem::UserType + 1 && items()[i]->type() == QGraphicsItem::UserType + 2)
          inter = static_cast<qPolyCurve*>(items()[i])
                      ->pointsOfIntersection(*static_cast<Bezier::Curve*>(static_cast<qCurve*>(items()[k])));
        if (items()[k]->type() == QGraphicsItem::UserType + 2 && items()[i]->type() == QGraphicsItem::UserType + 1)
          inter = static_cast<qPolyCurve*>(items()[k])
                      ->pointsOfIntersection(*static_cast<Bezier::Curve*>(static_cast<qCurve*>(items()[i])));
        if (items()[k]->type() == QGraphicsItem::UserType + 2 && items()[i]->type() == QGraphicsItem::UserType + 2)
          inter = static_cast<qPolyCurve*>(items()[i])
                      ->pointsOfIntersection(*static_cast<Bezier::PolyCurve*>(static_cast<qPolyCurve*>(items()[k])));

        for (auto& dot : inter)
          painter->drawEllipse(QPointF(dot.x(), dot.y()), 3, 3);
      }
  }
}

void CustomScene::mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  const int sensitivity = 5;
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
  if (mouseEvent->button() == Qt::RightButton)
  {
    dot = addEllipse(QRectF(QPointF(p.x(), p.y()), QSizeF(6, 6)), QPen(Qt::yellow), QBrush(Qt::red, Qt::SolidPattern));
    for (auto&& curve : items())
    {
      if (is_curve)
      {
        auto t1 = c_curve->projectPoint(p);
        auto p1 = c_curve->valueAt(t1);
        auto tan1 = c_curve->tangentAt(t1);
        line.insert(curve, addLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())), QPen(Qt::red)));
        tan.insert(curve, addLine(QLineF(QPointF(p1.x(), p1.y()) - 150 * QPointF(tan1.x(), tan1.y()),
                                         QPointF(p1.x(), p1.y()) + 150 * QPointF(tan1.x(), tan1.y())),
                                  QPen(Qt::blue)));
        auto t2 = c_curve->iterateByLength(t1, 50);
        auto a = c_curve->valueAt(t2);
        byLength.insert(curve, addEllipse(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)), QPen(Qt::yellow), QBrush(Qt::red, Qt::SolidPattern)));
      }
      if (is_poly)
      {
        auto t1 = c_poly->projectPoint(p);
        auto p1 = c_poly->valueAt(t1);
        auto tan1 = c_poly->tangentAt(t1);
        line.insert(curve, addLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())), QPen(Qt::red)));
        tan.insert(curve, addLine(QLineF(QPointF(p1.x(), p1.y()) - 150 * QPointF(tan1.x(), tan1.y()),
                                         QPointF(p1.x(), p1.y()) + 150 * QPointF(tan1.x(), tan1.y())),
                                  QPen(Qt::blue)));
        auto t2 = c_poly->iterateByLength(t1, 50);
        auto a = c_poly->valueAt(t2);
        byLength.insert(curve, addEllipse(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)), QPen(Qt::yellow), QBrush(Qt::red, Qt::SolidPattern)));
      }
    }
    show_projection = true;
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    for (auto&& curve : items())
    {
      if (!is_curve && !is_poly)
        continue;

      if (mouseEvent->modifiers().testFlag(Qt::ControlModifier))
      {
        if (is_curve)
        {
          double t = c_curve->projectPoint(p);
          auto pt = c_curve->valueAt(t);
          if ((pt - p).norm() < 10)
            curve->setSelected(true);
        }
        if (is_poly)
        {
          double t = c_poly->projectPoint(p);
          auto pt = c_poly->valueAt(t);
          if ((pt - p).norm() < 10)
            curve->setSelected(true);
        }
      }
      else
      {
        for (auto&& item : selectedItems())
          item->setSelected(false);
        if (is_curve)
        {
          auto pv = c_curve->controlPoints();
          for (unsigned k = 0; k < pv.size(); k++)
            if ((pv[k] - p).norm() < sensitivity && c_curve->getDraw_control_points())
            {
              update_cp = true;
              cp_to_update = std::make_pair(curve, k);
            }
        }
        if (is_poly)
        {
          auto pv = c_poly->controlPoints();
          for (unsigned k = 0; k < pv.size(); k++)
            if ((pv[k] - p).norm() < sensitivity && c_poly->getDraw_control_points())
            {
              update_cp = true;
              cp_to_update = std::make_pair(curve, k);
            }
        }
        if (update_cp)
          break;
        if (is_curve)
        {
          double t = c_curve->projectPoint(p);
          auto pt = c_curve->valueAt(t);
          auto ep = c_curve->endPoints();
          if ((pt - p).norm() < 10 && (pt - ep.first).norm() > 20 && (pt - ep.second).norm() > 20)
          {
            update_curvature = true;
            t_to_update = std::make_pair(c_curve, t);
            break;
          }
        }
      }
    }
  }
  if (mouseEvent->button() == Qt::MiddleButton)
  {
    for (auto&& curve : items())
    {
      if (is_curve)
      {
        auto t = c_curve->projectPoint(p);
        if ((c_curve->valueAt(t) - p).norm() < sensitivity)
        {
          this->removeItem(curve);
          auto split = c_curve->splitCurve(t);
          delete curve;
          qCurve *c1, *c2;
          c1 = new qCurve(split.first);
          c2 = new qCurve(split.second);
          this->addItem(c1);
          this->addItem(c2);
          update();
          break;
        }
      }
    }
  }
}

void CustomScene::mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
  if (show_projection)
  {
    dot->setRect(QRectF(QPointF(p.x() - 3, p.y() - 3), QSizeF(6, 6)));
    for (auto&& curve : items())
    {
      if (is_curve)
      {
        auto t1 = c_curve->projectPoint(p);
        auto p_c = c_curve->valueAt(t1);
        auto p1 = c_curve->valueAt(t1);
        auto tan1 = c_curve->tangentAt(t1);
        line[curve]->setLine(QLineF(QPointF(p.x(), p.y()), QPointF(p_c.x(), p_c.y())));
        tan[curve]->setLine(QLineF(QPointF(p1.x(), p1.y()) - 500 * QPointF(tan1.x(), tan1.y()),
                                   QPointF(p1.x(), p1.y()) + 500 * QPointF(tan1.x(), tan1.y())));
        auto t2 = c_curve->iterateByLength(t1, 50);
        auto a = c_curve->valueAt(t2);
        byLength[curve]->setRect(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)));
      }
      else if (is_poly)
      {
        auto t1 = c_poly->projectPoint(p);
        auto p1 = c_poly->valueAt(t1);
        auto tan1 = c_poly->tangentAt(t1);
        line[curve]->setLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())));
        tan[curve]->setLine(QLineF(QPointF(p1.x(), p1.y()) - 500 * QPointF(tan1.x(), tan1.y()),
                                   QPointF(p1.x(), p1.y()) + 500 * QPointF(tan1.x(), tan1.y())));
        auto t2 = c_poly->iterateByLength(t1, 50);
        auto a = c_poly->valueAt(t2);
        byLength[curve]->setRect(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)));
      }
    }
  }
  if (update_cp)
  {
    auto curve = cp_to_update.first;
    if (is_curve)
    {
      c_curve->prepareGeometryChange();
      c_curve->manipulateControlPoint(cp_to_update.second, p);
    }
    if (is_poly)
    {
      c_poly->prepareGeometryChange();
      c_poly->manipulateControlPoint(cp_to_update.second, p);
    }
    update();
  }
  if (update_curvature)
  {
    try
    {
      t_to_update.first->prepareGeometryChange();
      t_to_update.first->manipulateCurvature(t_to_update.second, p);
      update();
    }
    catch (char const* err)
    {
      update_curvature = false;
      QMessageBox::warning(nullptr, "Warning", QString().sprintf("%s", err));
    }
  }
  QGraphicsScene::mouseMoveEvent(mouseEvent);
}

void CustomScene::mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  if (mouseEvent->button() == Qt::RightButton)
  {
    if (show_projection)
    {
      removeItem(dot);
      for (auto&& curve : items())
        if (is_curve || is_poly)
        {
          removeItem(line[curve]);
          removeItem(tan[curve]);
          removeItem(byLength[curve]);
        }
      line.clear();
      tan.clear();
      byLength.clear();
      show_projection = false;
    }
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    update_cp = false;
    update_curvature = false;
  }
}

void CustomScene::keyPressEvent(QKeyEvent* keyEvent)
{
  if (keyEvent->key() == 72) // key H
  {
    QMessageBox::information(nullptr, "Help", "\
Mouse controls:\n\
Right click - project mouse pointer on all curves\n\
Ctrl + Scroll - zoom in/out\n\
Ctrl + Left click - select curves\n\
Left click - deselect curves\n\
Left click + drag control point - manipulate control point\n\
Left click + drag curve - manipulate curve (only for 2nd and 3rd order)\n\
\n\
Keyboard shortcuts:\n\
H - display help\n\
B - toggle bounding box display\n\
I - toggle intesections display\n\
C - toggle curvature display (of selected curves)\n\
P - toggle control points display (of selected curves)\n\
Key Up - raise the order (of selected curves)\n\
Key Down - lower the order (of selected curves)\n\
Key + - join multiple curves into polycurve\n\
Delete - delete curve/polycurve");
  }

  if (keyEvent->key() == 66) // key B
  {
    draw_box_ = !draw_box_;
    update();
  }
  if (keyEvent->key() == 73) // key I
  {
    draw_inter_ = !draw_inter_;
    update();
  }
  if (keyEvent->key() == 67) // key C
  {
    for (auto&& curve : selectedItems())
      if (is_curve)
        c_curve->setDraw_curvature_radious(!c_curve->getDraw_curvature_radious());
      else if (is_poly)
        c_poly->setDraw_curvature_radious(!c_poly->getDraw_curvature_radious());
    update();
  }
  if (keyEvent->key() == 80) // key P
  {
    for (auto&& curve : selectedItems())
      if (is_curve)
        c_curve->setDraw_control_points(!c_curve->getDraw_control_points());
      else if (is_poly)
        c_poly->setDraw_control_points(!c_poly->getDraw_control_points());
    update();
  }
  if (keyEvent->key() == 16777235) // key UP
  {
    for (auto&& curve : selectedItems())
      if (is_curve)
        c_curve->elevateOrder();
    update();
  }
  if (keyEvent->key() == 16777237) // key DOWN
  {
    for (auto&& curve : selectedItems())
      if (is_curve)
        try
        {
          c_curve->lowerOrder();
        }
        catch (char const* err)
        {
          QMessageBox::warning(nullptr, "Warning", QString().sprintf("%s", err));
        }
    update();
  }
  if (keyEvent->key() == 43) // key +
  {
    qPolyCurve* new_poly = nullptr;
    for (auto&& curve : selectedItems())
      if (is_curve)
      {
        std::shared_ptr<Bezier::Curve> ptr(c_curve);
        if (new_poly)
          new_poly->insertBack(ptr);
        else
          new_poly = new qPolyCurve(ptr);
        removeItem(curve);
      }
    if (new_poly)
      addItem(new_poly);
    update();
  }
  if (keyEvent->key() == 16777223) // Delete
  {
    for (auto&& curve : selectedItems())
    {
      removeItem(curve);
    }
    update();
  }
  if (keyEvent->key() == Qt::Key_L) {
      for (auto&& curve : selectedItems()) {
          if (is_curve) {
              c_curve->setLocked(!c_curve->getLocked());
          }
      }
      update();
  }
  if (keyEvent->key() == Qt::Key_Enter || keyEvent->key() == Qt::Key_Return) {
      if (selectedItems().size() != 2) {
          return;
      }

      qCurve* curve_locked = nullptr;
      qCurve* curve_unlocked = nullptr;

      if (c_curve_item(selectedItems().first())->getLocked()) {
          curve_locked = c_curve_item(selectedItems().first());
          curve_unlocked = c_curve_item(selectedItems().last());
      }
      else if (c_curve_item(selectedItems().last())->getLocked()) {
          curve_locked = c_curve_item(selectedItems().last());
          curve_unlocked = c_curve_item(selectedItems().first());
      }

      if (curve_locked != nullptr && curve_unlocked != nullptr) {
          while (curve_locked->order() < 5) {
              curve_locked->elevateOrder();
          }
          while (curve_unlocked->order() < 5) {
              curve_unlocked->elevateOrder();
          }

          std::vector<double> beta_coeffs = {1, 0, 0};
          curve_unlocked->applyContinuity(*dynamic_cast<Bezier::Curve*>(curve_locked), beta_coeffs);
          update();
      }
  }

  //  if (keyEvent->key() >= 48 && keyEvent->key() <= 59) // num keys
  //  {
  //    Bezier::Continuity c;
  //    c.order = keyEvent->key() - 48;
  //    c.type = keyEvent->modifiers().testFlag(Qt::ControlModifier) ? 'G' : 'C';
  //    for (auto&& curve : selectedItems())
  //    {
  //      if(is_poly)
  //        for(unsigned k = 0; k < c_poly->getSize() - 1; k++)
  //          c_poly->setContinuity(k, c);
  //    }
  //    update();
  //  }
}
